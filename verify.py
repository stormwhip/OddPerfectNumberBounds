"""Verify q > 220000000000000000000 for restrained odd perfect numbers.
Run: python verify.py. Requires sieve.tsv in the same directory.
All certificate comparisons use exact rational arithmetic.
"""
import csv
from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor, prod
from pathlib import Path
S=10**18
P=(7,13,19,31,37,43,61,67,73,79,97,103)
T=[(7,a,b) for a,b in combinations(P[1:],2) if 7*a*b<=10000]
B=50290000000000
Z=57623000000000000000000
Q=220000000000000000000
def Log(x):
    t=(x-1)/(x+1)
    y=t*t
    p=t
    l=F(0)
    for i in range(24):
        l+=2*p/(2*i+1)
        p*=y
    return l,l+2*p/(49*(1-y))
def Bounds(x):
    x=F(x)
    n=0
    while x>2:
        x/=2
        n+=1
    l,r=Log(x)
    a,b=Log(F(2))
    return F(floor((l+n*a)*S),S),F(ceil((r+n*b)*S),S)
def Certificate():
    a=[(c,(-1)**k) for k in range(3) for c in combinations(P,k)]
    a.extend((c,-1) for c in T)
    assert len(a)==97 and len(T)==18
    assert all(prod(c)<=10000 for c,_ in a)
    for mask in range(1<<len(P)):
        b={p:bool(mask&(1<<i)) for i,p in enumerate(P)}
        k=mask.bit_count()
        g=1-k+k*(k-1)//2
        g-=sum(all(b[p] for p in c) for c in T)
        assert g>=(1 if mask==0 else 0)
    d=sum((s*prod(F(2,p-1) for p in c) for c,s in a),F(0))
    e=sum((F(2**len(c),840) for c,_ in a),F(0))
    assert d==F(64959667,169646400) and e==F(433,840)
    return d,e
def Sieve():
    with Path(__file__).with_name('sieve.tsv').open(newline='',encoding='utf-8') as f:
        a=list(csv.DictReader(f,delimiter='\t'))
    assert a and int(a[0]['lower_exclusive'])==3000000
    n=s=0
    l=3000000
    for x in a:
        assert int(x['lower_exclusive'])==l
        l=int(x['upper_inclusive'])
        n+=int(x['survivors'])
        s+=int(x['scaled_reciprocal_upper'])
    assert l==B
    return n,s,len(a)
def Verify(n,s,d,e):
    b,_=Bounds(B)
    _,z=Bounds(Z)
    _,v=Bounds(z/b)
    t=F(B,B-1)
    li=F(Z,1)/z-F(B,1)/b+F(Z-B,1)/z**2
    m=d*(t*v-li/(Z-1))
    err=e*(t/b**2+t*t*(1/b-1/z))
    ans=F(s,S)+m+err+F(2*Q+2-n,Z-1)+F(1,2999999)
    target,_=Bounds(F(4,3))
    return ans,target
def Solve():
    d,e=Certificate()
    e-=F(1,840)-F(8375,10000000)
    e-=2*(F(1,840)-F(486,1000000))
    assert e==F(21575999,42000000)
    n,s,c=Sieve()
    assert (n,s,c)==(189326366793,86154121237231668,4943)
    ans,target=Verify(n,s,d,e)
    scale=10**15
    u=ceil(ans*scale)
    v=floor(target*scale)
    print(f'Validated {1<<len(P)} Boolean assignments and {c} sieve rows.')
    print(f'97 majorant terms; {len(T)} selected triple intersections.')
    print(f'delta = {d}; C = {e}')
    print(f'upper <= {u//scale}.{u%scale:015d}')
    print(f'log(4/3) >= {v//scale}.{v%scale:015d}')
    if ans>=target:
        raise ArithmeticError('Final inequality does not certify the bound.')
    print(f'CERTIFIED: q > {Q}')
if __name__=='__main__':
    Solve()
