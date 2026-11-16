"""Deterministic constraint-generation baseline for Stage 2 numerical LPs.

Run from repository root:
python "Stage 2/experiments/constraint_generation.py" --primes 15
All Boolean states are separated explicitly at the end of each iteration;
this is a floating-point experiment, not a number-theory certificate.
"""
import argparse, json, time
from itertools import combinations
from math import prod, log
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack
from extended_lp_search import PRIMES, B

Z = 57433500000000000000000

def build(primes, degree, maxmod):
    n = len(primes)
    terms = [(sub, prod(primes[i] for i in sub)) for k in range(degree+1)
             for sub in combinations(range(n), k)
             if prod(primes[i] for i in sub) <= maxmod]
    masks = np.arange(1, 1<<n, dtype=np.uint32)
    rr, cc = [], []
    for j, (sub, _) in enumerate(terms):
        bits = sum(1<<i for i in sub)
        rows = np.flatnonzero((masks & bits) == bits)
        rr.append(rows)
        cc.append(np.full(len(rows),j,dtype=np.int32))
    rows, cols = np.concatenate(rr), np.concatenate(cc)
    M = csr_matrix((np.ones(len(rows)),(rows,cols)),
                   shape=(len(masks), len(terms)))
    den = np.array([prod(2/(primes[i]-1) for i in sub) for sub,_ in terms])
    err = np.array([2**len(sub)/(840 if mod<=10000 else 160)
                    for sub,mod in terms])
    return M, terms, den, err

def solve(n=15, degree=4, maxmod=100000, batch=128, maxrounds=80):
    ps = PRIMES[:n]
    start = time.perf_counter()
    M, terms, den, err = build(ps, degree, maxmod)
    b,z = log(B),log(Z)
    rho = B/(B-1)
    A = rho*log(z/b)-(Z/z-B/b+(Z-B)/z**2)/(Z-1)
    D = rho/b**2+rho**2*(1/b-1/z)
    c = np.r_[A*den+D*err,-A*den+D*err]
    l = len(terms)
    states = np.arange(1,1<<n,dtype=np.uint32)
    selected = np.array(sorted(set(sum(1<<i for i in sub)-1
                                   for sub,_ in terms if sub)),dtype=np.int64)
    added = np.zeros(len(states),dtype=bool)
    added[selected] = True
    history = []
    for step in range(maxrounds):
        W = M[selected,:]
        fit = linprog(c,A_ub=hstack((-W,W),format='csr'),
                      b_ub=np.zeros(len(selected)),
                      bounds=[(1,1)]+[(0,None)]*(l-1)+[(0,0)]+[(0,None)]*(l-1),
                      method='highs',options={'time_limit':25})
        if not fit.success:
            return {'status':fit.message,'history':history}
        a = fit.x[:l]-fit.x[l:]
        values = M@a
        violated = np.flatnonzero((values < -1e-8)&(~added))
        history.append({'step':step,'active':len(selected),
                        'violations':len(violated),
                        'worst':float(np.min(values))})
        if len(violated)==0:
            return {'status':'no_boolean_violations','primes':list(ps),
                    'columns':l,'states':len(states),
                    'active_constraints':len(selected),'rounds':step+1,
                    'seconds':round(time.perf_counter()-start,4),
                    'delta':float(a@den),'C':float(abs(a)@err),
                    'objective':float(fit.fun),
                    'nonzero':int(np.sum(abs(a)>1e-7)),
                    'history':history,'certified':False}
        worst = violated[np.argsort(values[violated])[:batch]]
        selected = np.unique(np.r_[selected,worst])
        added[worst] = True
    return {'status':'round_limit','history':history}

if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--primes',type=int,choices=range(6,17),default=15)
    p.add_argument('--batch',type=int,default=128)
    p.add_argument('--output')
    a = p.parse_args()
    result = solve(n=a.primes,batch=a.batch)
    print(json.dumps({k:v for k,v in result.items() if k!='history'},indent=2))
    if a.output:
        from pathlib import Path
        file = Path(a.output)
        file.parent.mkdir(parents=True,exist_ok=True)
        file.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
