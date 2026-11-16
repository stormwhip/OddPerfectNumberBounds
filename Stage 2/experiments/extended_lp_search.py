"""Explore degree <=4 inclusion-exclusion LPs for a fixed B and Z.

Local numerical optimization only; a new number-theory theorem requires
independent exact rational verification. No API or paid models.
Requires numpy and scipy. Run from the repository root:
python "Stage 2/experiments/extended_lp_search.py" --primes 14 --degree 4 --maxmod 100000
"""
import argparse
import json
from itertools import combinations
from math import prod, log
from time import perf_counter
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack

PRIMES=(7,13,19,31,37,43,61,67,73,79,97,103,109,127,139,151)
B=50290000000000
DEFAULT_Z=57433500000000000000000


def solve(primes, degree, maxmod, Z):
    n=len(primes)
    terms=[(subset,prod(primes[j] for j in subset))
           for k in range(degree+1)
           for subset in combinations(range(n),k)
           if prod(primes[j] for j in subset)<=maxmod]
    masks=np.arange(1,1<<n,dtype=np.uint32)
    rr=[]; cc=[]
    for j,(subset,_) in enumerate(terms):
        bits=sum(1<<i for i in subset)
        rows=np.flatnonzero((masks&bits)==bits)
        rr.append(rows)
        cc.append(np.full(len(rows),j,dtype=np.int32))
    rows=np.concatenate(rr); cols=np.concatenate(cc)
    M=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(masks),len(terms)))
    den=np.array([prod(2/(primes[j]-1) for j in subset) for subset,_ in terms])
    err=np.array([2**len(subset)/(840 if m<=10000 else 160)
                  for subset,m in terms])
    b,z=log(B),log(Z)
    rho=B/(B-1)
    A=rho*log(z/b)-(Z/z-B/b+(Z-B)/z**2)/(Z-1)
    D=rho/b**2+rho**2*(1/b-1/z)
    costs=np.r_[A*den+D*err,-A*den+D*err]
    m=len(terms)
    start=perf_counter()
    result=linprog(costs,A_ub=hstack((-M,M),format='csr'),
                   b_ub=np.zeros(M.shape[0]),
                   bounds=[(1,1)]+[(0,None)]*(m-1)+[(0,0)]+[(0,None)]*(m-1),
                   method='highs',options={'time_limit':35})
    output=dict(primes=list(primes),degree=degree,maxmod=maxmod,Z=str(Z),
                columns=m,boolean_constraints=M.shape[0],matrix_nonzeros=M.nnz,
                seconds=round(perf_counter()-start,4),status=result.message,
                numerical_only=True,certified=False)
    if result.x is not None:
        coef=result.x[:m]-result.x[m:]
        delta=float(coef@den); C=float(abs(coef)@err)
        output.update(delta=delta,C=C,objective=A*delta+D*C,
                      nonzero=int(sum(abs(coef)>1e-7)),
                      higher_order_terms=[dict(primes=[primes[i] for i in sub],
                                              coefficient=round(float(coef[j]),7))
                        for j,(sub,_) in enumerate(terms)
                        if len(sub)>=3 and abs(coef[j])>1e-7])
    return output


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--primes',type=int,choices=range(6,17),default=14)
    ap.add_argument('--degree',type=int,choices=range(1,5),default=4)
    ap.add_argument('--maxmod',type=int,choices=(10000,100000),default=100000)
    ap.add_argument('--Z',type=int,default=DEFAULT_Z)
    ap.add_argument('--output',default=None)
    args=ap.parse_args()
    if args.Z <= B:ap.error('Z must exceed B')
    result=solve(PRIMES[:args.primes],args.degree,args.maxmod,args.Z)
    raw=json.dumps(result,ensure_ascii=False,indent=2)
    print(raw)
    if args.output:
        from pathlib import Path
        p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(raw+'\n',encoding='utf8')


if __name__=='__main__':main()
