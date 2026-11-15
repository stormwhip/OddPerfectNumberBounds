"""Run complete sieve intervals until a deadline, then verify and update the report."""

import argparse
import csv
import ctypes
import io
import os
import re
import subprocess
import tempfile
import time
import traceback
from datetime import datetime
from itertools import groupby
from math import ceil, floor
from pathlib import Path

from QBoundVerify import verify


SEED_LIMIT = 3_000_000
UPPER_LIMIT = 10**14  # Must match the validated range in QBoundPrefix.cpp.
SCALE = 10**18
HEADERS = (
    "lower_exclusive", "upper_inclusive", "residue_limit",
    "survivors", "scaled_reciprocal_upper",
)


def now():
    return datetime.now().astimezone()


def log(message):
    print(f"{now().isoformat(timespec='seconds')} {message}", flush=True)


def write_atomic(path, content):
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as output:
        temporary = Path(output.name)
        output.write(content)
        output.flush()
        os.fsync(output.fileno())
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def read_rows(path):
    with path.open(encoding="utf-8", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        if tuple(next(reader)) != HEADERS:
            raise ValueError("Unexpected sieve data header")
        rows = [tuple(map(int, row)) for row in reader]
    previous = SEED_LIMIT
    for lower, upper, residue_limit, count, scaled in rows:
        if not (lower == previous < upper <= UPPER_LIMIT
                and 0 <= residue_limit <= SEED_LIMIT and 0 <= count <= upper - lower
                and count * SCALE <= scaled * (upper - 1)
                and scaled <= count * ((SCALE + lower - 1) // lower)):
            raise ValueError(f"Invalid interval or result: {lower}-{upper}")
        previous = upper
    if not rows:
        raise ValueError("No prior sieve checkpoints")
    return rows


def save_rows(path, rows):
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADERS)
    writer.writerows(rows)
    write_atomic(path, output.getvalue())


def run_sieve(args, path, rows):
    while rows[-1][1] < UPPER_LIMIT:
        seconds = (args.until - now()).total_seconds()
        if seconds <= 0:
            return
        lower = rows[-1][1]
        upper = min(lower + args.step, UPPER_LIMIT)
        log(f"Sifting ({lower}, {upper}]")
        try:
            output = subprocess.run(
                [str(args.executable), str(upper), str(args.threads),
                 str(lower), str(SEED_LIMIT)],
                check=True, capture_output=True, text=True, timeout=seconds,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
        except subprocess.TimeoutExpired:
            log("Deadline reached during an interval")
            return
        if now() >= args.until:
            log("Deadline reached before checkpoint")
            return
        count, scaled = map(int, output.stdout.split())
        if not (0 <= count <= upper - lower and count * SCALE <= scaled * (upper - 1)
                and scaled <= count * ((SCALE + lower - 1) // lower)):
            raise ValueError(f"Unexpected sieve output for endpoint {upper}")
        rows.append((lower, upper, SEED_LIMIT, count, scaled))
        save_rows(path, rows)
        log(f"Saved B={upper}: {count} candidates, scaled sum {scaled}")
    log(f"Reached validated arithmetic limit B={UPPER_LIMIT}")


def certified_bound(B, count, scaled):
    def capacity(Z):
        total, target = verify(B, Z, 0, count, scaled)
        return ceil((target - total) * (Z - 1) / 2) - 1

    _, coarse_z = max(
        (capacity(Z), Z)
        for decade in range(6, 11) for digit in range(1, 10)
        for Z in (B * digit * 10**decade,)
    )
    maximum, Z = max((capacity(z), z) for percent in range(50, 151)
                     for z in (coarse_z * percent // 100,))
    if maximum <= 0:
        raise ValueError("Could not certify a positive bound")
    unit = 10 ** max(0, len(str(maximum)) - 4)
    q = maximum // unit * unit
    total, target = verify(B, Z, q, count, scaled)
    if total >= target:
        raise ValueError("Numerical inequality failed")
    return q, Z, total, target


def spaced(number):
    return f"{number:,}".replace(",", " ")


def latex_integer(number):
    return f"{number:,}".replace(",", r"\,")



def update_report(path, rows, started, deadline, failure):
    report = path.read_text(encoding="utf-8")
    B = rows[-1][1]
    count = sum(row[3] for row in rows)
    scaled = sum(row[4] for row in rows)
    q, Z, total, target = certified_bound(B, count, scaled)
    digits = 10**12
    upper, lower = ceil(total * digits), floor(target * digits)
    upper_text = f"{upper // digits}.{upper % digits:012d}"
    lower_text = f"{lower // digits}.{lower % digits:012d}"
    q_digits = str(q)
    coefficient = (q_digits[0] + "." + q_digits[1:4]).rstrip("0").rstrip(".")
    q_latex = rf"{coefficient}\times10^{{{len(q_digits) - 1}}}"
    ratio_coefficient, ratio_exponent = f"{q / 23_000_000_000:.2e}".split("e")

    def field_pattern(name):
        return re.compile(rf"(<!--{name}-->)(.*?)(<!--/{name}-->)", re.S)

    def get_field(name):
        match = field_pattern(name).search(report)
        if match is None:
            raise ValueError(f"Missing report field: {name}")
        return match[2]

    def set_field(name, value):
        nonlocal report
        report, matches = field_pattern(name).subn(
            lambda match: match[1] + value + match[3], report
        )
        if matches != 1:
            raise ValueError(f"Expected one report field: {name}")

    # Keep localized wording in the report and update only marked data.
    date_parts = iter(map(str, (deadline.year, deadline.month, deadline.day)))
    set_field("ReportDate", re.sub(r"\d+", lambda _: next(date_parts),
                                  get_field("ReportDate"), count=3))
    set_field("SieveEndpoint", rf"$B={latex_integer(B)}$")
    set_field("Result", f"\n\n$$\n\\boxed{{q>{q_latex}}}.\n$$\n\n")
    set_field("Improvement", rf"${ratio_coefficient}\times10^{{{int(ratio_exponent)}}}$")
    set_field("Started", started.strftime("%Y-%m-%d %H:%M"))
    set_field("Deadline", deadline.strftime("%Y-%m-%d %H:%M"))
    set_field("RunNote", " [Run log](QBoundOvernight.log) " if failure else "")
    set_field("Parameters", rf"$B={latex_integer(B)}, Z={latex_integer(Z)}$")
    set_field("ScaledSum", rf"$S={scaled // SCALE}.{scaled % SCALE:018d}$")
    set_field("CandidateCount", rf"$K_B={latex_integer(count)}$")
    set_field("TrialBound", rf"$q\le {latex_integer(q)}$")

    previous_table = get_field("SummaryTable").strip().splitlines()
    table = previous_table[:2]
    total_label = previous_table[-1].split("|")[1].strip()
    for limit, group in groupby(rows, key=lambda row: row[2]):
        group = list(group)
        lo, hi = group[0][0], group[-1][1]
        n = sum(row[3] for row in group)
        s = sum(row[4] for row in group)
        interval = rf"$({latex_integer(lo)},{latex_integer(hi)}]$"
        table.append(f"| {interval} | {spaced(limit)} | {spaced(n)} | {spaced(s)} |")
    table.append(f"| {total_label} |  | {spaced(count)} | {spaced(scaled)} |")
    set_field("SummaryTable", "\n\n" + "\n".join(table) + "\n\n")

    report, upper_matches = re.subn(r"(?m)^(\\le )\d+\.\d+\.$",
                                    lambda match: match[1] + upper_text + ".", report)
    report, lower_matches = re.subn(r"(?<=\\log\(4/3\)\\ge )\d+\.\d+", lower_text, report)
    if (upper_matches, lower_matches) != (1, 1):
        raise ValueError("Could not locate the numerical comparison")
    commands = iter((f"{count} {scaled}", "$count $scaledSum"))
    report, matches = re.subn(
        r"(?m)^python QBoundVerify\.py .+$",
        lambda _: f"python QBoundVerify.py {B} {Z} {q} {next(commands)}", report,
    )
    if matches != 2:
        raise ValueError("Could not locate both verification commands")
    write_atomic(path, report)
    log(f"Report updated: B={B}, q>{q}, upper={upper_text}, target>={lower_text}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--until", type=datetime.fromisoformat, required=True)
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--step", type=int, default=10**10)
    args = parser.parse_args()
    started = now()
    if (args.until.tzinfo is None or args.until <= started
            or not 1 <= args.threads <= 128 or args.step <= 0):
        parser.error("supply a future deadline with UTC offset and valid thread and step counts")
    args.executable = args.executable.resolve(strict=True)
    data = args.workspace / "QBoundSieveRuns.tsv"
    report = args.workspace / "ResearchReport.md"
    rows = read_rows(data)
    execution_state = None
    if os.name == "nt":
        execution_state = ctypes.windll.kernel32.SetThreadExecutionState
        execution_state.argtypes = [ctypes.c_uint32]
        execution_state.restype = ctypes.c_uint32
        if not execution_state(0x80000001):
            log("Windows could not hold the system awake")
    log(f"Started at B={rows[-1][1]}, deadline={args.until.isoformat()}, threads={args.threads}")
    failure = False
    try:
        try:
            run_sieve(args, data, rows)
        except (OSError, subprocess.SubprocessError, ValueError):
            traceback.print_exc()
            failure = True
        while (remaining := (args.until - now()).total_seconds()) > 0:
            time.sleep(min(remaining, 60))
        # Reload only the committed, contiguous checkpoint data.
        update_report(report, read_rows(data), started, args.until, failure)
    finally:
        if execution_state is not None:
            execution_state(0x80000000)


if __name__ == "__main__":
    main()
