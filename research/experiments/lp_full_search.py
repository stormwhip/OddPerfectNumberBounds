"""Explore the finite Boolean LP for fixed Z, without ML or paid services.

Run: python research/experiments/lp_full_search.py 12 10000
     python research/experiments/lp_full_search.py 12 100000
     python research/experiments/lp_full_search.py 13 100000
Floating-point LP optima are exploratory. Certify any new bound separately.
"""
import argparse
import time
from itertools import combinations
from math import prod, log
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack

POOL=(7,13,19,31,37,43,61,67,73,79,97,103,109,127)
B,Z=50290000000000,57433500000000000000000


def run(n, maxmod):
    P=POOL[:n]
    b,z=log(B),log(Z)
    rho=B/(B-1)
    A=rho*log(z/b)-(Z/z-B/b+(Z-B)/z**2)/(Z-1)
    D=rho/b**2+rho*rho*(1/b-1/z)
    terms=[(sub,prod(P[j] for j in sub))
           for k in range(4) for sub in combinations(range(n),k)
           if prod(P[j] for j in sub)<=maxmod]
    m=len(terms)
    masks=np.arange(1,1<<n,dtype=np.uint16)
    rr,cc=[],[]
    for j,(sub,_) in enumerate(terms):
        bits=sum(1<<i for i in sub)
        rows=np.flatnonzero((masks&bits)==bits)
        rr.extend(rows.tolist())
        cc.extend([j]*len(rows))
    M=csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(masks),m))
    objective_matrix=hstack((-M,M),format="csr")
    density=np.array([prod(2/(P[j]-1) for j in sub) for sub,_ in terms])
    errors=np.array([2**len(sub)/(840 if modulus<=10000 else 160)
                     for sub,modulus in terms])
    objective=np.r_[A*density+D*errors,-A*density+D*errors]
    t=time.monotonic()
    result=linprog(objective,A_ub=objective_matrix,b_ub=np.zeros(len(masks)),
                   bounds=[(1,1)]+[(0,None)]*(m-1)+[(0,0)]+[(0,None)]*(m-1),
                   method="highs",options={"time_limit":20})
    assert result.success,result.message
    coef=result.x[:m]-result.x[m:]
    print(f"n={n} maxmod={maxmod} columns={m} "
          f"seconds={time.monotonic()-t:.3f} status=optimal")
    print(f"delta={coef@density:.15f} C={np.abs(coef)@errors:.15f} "
          f"nonzero={sum(np.abs(coef)>1e-7)} "
          f"objective={A*(coef@density)+D*(np.abs(coef)@errors):.15f}")
    print("nonzero cubic monomials:",
          [(tuple(P[j] for j in sub),round(float(coef[i]),5))
           for i,(sub,mod) in enumerate(terms)
           if len(sub)==3 and abs(coef[i])>1e-7])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("n",type=int,choices=range(6,15))
    p.add_argument("maxmod",type=int,choices=(10000,100000))
    a=p.parse_args()
    run(a.n,a.maxmod)
