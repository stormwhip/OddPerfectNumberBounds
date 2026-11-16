"""Local LightGBM vs LP pilot for certified sieve-majorant column selection.

Dependencies: numpy, scipy, scikit-learn, lightgbm.
No network, credentials, or paid APIs are used. Floating-point LP values
are experimental objectives, NOT new certified number-theory bounds.
"""
import argparse
import json
import warnings
from itertools import combinations
from math import prod, log
from pathlib import Path

import numpy as np
from lightgbm import LGBMClassifier
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore", message="X does not have valid feature names")
POOL = (7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103)
B, Z = 50290000000000, 57433500000000000000000
b, z = log(B), log(Z)
rho = B/(B-1)
A = rho*log(z/b)-(Z/z-B/b+(Z-B)/z**2)/(Z-1)
D = rho/b**2+rho*rho*(1/b-1/z)


def features(triple, primes):
    t = prod(triple)
    return [int(7 in triple), log(t), sum(1/(p-1) for p in triple),
            prod(2/(p-1) for p in triple), min(triple), max(triple),
            sum(triple), len(primes), sum(p <= 43 for p in primes)]


def instance(primes):
    n = len(primes)
    base = [()] + list(combinations(range(n), 1)) + list(combinations(range(n), 2))
    extras = [sub for sub in combinations(range(n), 3)
              if prod(primes[i] for i in sub) <= 10000]
    terms = base+extras
    masks = np.arange(1, 1 << n, dtype=np.uint16)
    rows, cols = [], []
    for j, sub in enumerate(terms):
        mask = sum(1 << i for i in sub)
        nz = np.flatnonzero((masks & mask) == mask)
        rows.extend(nz.tolist())
        cols.extend([j]*len(nz))
    mat = csr_matrix((np.ones(len(rows)), (rows, cols)),
                     shape=(len(masks), len(terms)))
    density = np.array([prod(2/(primes[i]-1) for i in sub) for sub in terms])
    errors = np.array([2**len(sub)/840 for sub in terms])
    feats = np.array([features(tuple(primes[i] for i in sub), primes)
                      for sub in extras]).reshape(-1, 9)
    return mat, density, errors, len(base), extras, feats


def solve(mat, density, errors, columns):
    ids = [0]+list(columns)
    m = len(ids)
    d, e = density[ids], errors[ids]
    objective = np.r_[A*d+D*e, -A*d+D*e]
    M = mat[:, ids]
    bounds = [(1, 1)] + [(0, None)]*(m-1) + [(0, 0)] + [(0, None)]*(m-1)
    out = linprog(objective, A_ub=hstack((-M, M), format="csr"),
                  b_ub=np.zeros(mat.shape[0]), bounds=bounds,
                  method="highs", options={"time_limit": 15})
    if not out.success:
        raise RuntimeError(out.message)
    a = out.x[:m]-out.x[m:]
    return float(A*(a@d)+D*(np.abs(a)@e)), a


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--train", type=int, default=60)
    cli.add_argument("--test", type=int, default=20)
    cli.add_argument("--budget", type=int, default=5)
    cli.add_argument("--seed", type=int, default=42)
    cli.add_argument("--output", default="Stage 2/experiments/pilot_results.json")
    args = cli.parse_args()
    rng = np.random.default_rng(args.seed)
    choices = [(7, *sub) for n in (8, 9, 10, 11)
               for sub in combinations(POOL[1:], n-1)]
    picked = rng.choice(len(choices), size=args.train+args.test, replace=False)
    cases, x_train, y_train = [], [], []
    for idx, choice in enumerate(picked):
        ps = choices[choice]
        mat, d, e, base, extras, x = instance(ps)
        full, coeffs = solve(mat, d, e, range(1, base+len(extras)))
        labels = (np.abs(coeffs[base:]) > 1e-7).astype(int)
        cases.append((ps, mat, d, e, base, extras, x, labels, full))
        if idx < args.train:
            x_train.extend(x.tolist())
            y_train.extend(labels.tolist())
    model = LGBMClassifier(n_estimators=70, num_leaves=9, max_depth=4,
                           min_child_samples=5, learning_rate=0.07,
                           random_state=args.seed, verbosity=-1, n_jobs=1)
    model.fit(np.asarray(x_train, dtype=float), np.asarray(y_train))
    records, all_labels, all_probs = [], [], []
    for ps, mat, d, e, base, extras, x, labels, full in cases[args.train:]:
        n = len(extras)
        probs = model.predict_proba(x)[:, 1] if n else np.array([])
        all_labels.extend(labels.tolist())
        all_probs.extend(probs.tolist())
        k = min(args.budget, n)
        selected_ml = np.argsort(-probs)[:k]
        heuristic = np.array([
            1000*(7 in tuple(ps[i] for i in sub))+
            A*d[base+j]-D*e[base+j] for j, sub in enumerate(extras)])
        selected_heur = np.argsort(-heuristic)[:k]
        selected_random = rng.choice(n, size=k, replace=False)
        base_columns = list(range(1, base))
        score0, _ = solve(mat, d, e, base_columns)
        def score(ids):
            return solve(mat, d, e, base_columns+[base+j for j in ids])[0]
        sm, sh, sr = score(selected_ml), score(selected_heur), score(selected_random)
        records.append({"prime_set": list(ps), "triples": n,
                        "gain_full": score0-full, "gain_ml": score0-sm,
                        "gain_heur": score0-sh, "gain_random": score0-sr,
                        "ml_selected_true": int(sum(labels[selected_ml]))})
    gain_names = ("gain_full", "gain_ml", "gain_heur", "gain_random")
    mean_gains = {name: float(np.mean([r[name] for r in records]))
                  for name in gain_names}
    ratios = {name: mean_gains[name]/mean_gains["gain_full"]
              for name in gain_names[1:]}
    output = {"seed": args.seed, "train_instances": args.train,
              "test_instances": args.test, "budget": args.budget,
              "model_average_precision": float(average_precision_score(all_labels, all_probs)),
              "model_roc_auc": float(roc_auc_score(all_labels, all_probs)),
              "average_tail_objective_gains": mean_gains,
              "fractions_of_full_lp_improvement": ratios,
              "ml_wins_over_rule": sum(r["gain_ml"]>r["gain_heur"]+1e-12
                                       for r in records),
              "cases": records, "proof_status": "exploratory float LP only"}
    dest = Path(args.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(output, indent=2, ensure_ascii=False)+"\n", encoding="utf8")
    print(json.dumps({key: val for key, val in output.items() if key != "cases"},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
