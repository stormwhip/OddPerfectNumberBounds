/*
Build: g++ -std=c++17 -O2 -pthread chain.cpp -o chain
Usage: ./chain [upper [threads [lower [second_limit]]]]
Output: two-step sieve counts, scaled reciprocal sums and sample witnesses.
*/
#include<bits/stdc++.h>
#define int long long
using namespace std;
using ull=unsigned long long;
using u128=__uint128_t;
const int N=1<<21;
const ull D=3000000,S=1000000000000000000ULL,MAXN=1000000000;
ull l=D,r=MAXN,t=1,lim=10000;

enum kind{base,back,forwardOrd,forwardSpec,twoSpec,cond,unresolved,forwardPrime,K};
const char* name[K]={"baseline","backward","forward_ordinary","forward_special",
    "two_specials","special_conditional","unresolved","forward_prime_quotient"};
struct root
{
	ull p,x,y;
};
struct witness
{
	ull r=0,p=0,a=0,b=0,q=0,c=0;
};
struct node
{
	ull cnt[K]{},sum[K]{},prime=0,composite=0;
	witness first[K];
	unordered_map<ull,ull> group;
};

ull Mul(ull a,ull b,ull mod)
{
	if(mod<=UINT32_MAX) return a*b%mod;
	return (u128)a*b%mod;
}
ull Pows(ull a,ull b,ull mod)
{
	ull ans=1;
	while(b)
	{
		if(b&1) ans=Mul(ans,a,mod);
		a=Mul(a,a,mod);
		b>>=1;
	}
	return ans;
}
bool IsPrime(ull n)
{
	if(n%2==0) return n==2;
	if(n<2) return false;
	ull d=n-1,s=0;
	while(!(d&1)) d>>=1,s++;
	// Sorenson and Webster's bound for the first twelve prime bases.
	for(ull a:{2,3,5,7,11,13,17,19,23,29,31,37})
	{
		if(a%n==0) continue;
		ull x=Pows(a%n,d,n),i=0;
		if(x==1) continue;
		for(;i<s&&x!=n-1;i++) x=Mul(x,x,n);
		if(i==s) return false;
	}
	return true;
}
ull Factor(ull n)
{
	for(ull c=1;c<=8;c++)
	{
		auto Next=[&](ull x){return (Mul(x,x,n)+c)%n;};
		ull x=0,y=2,ys=0,g=1;
		for(ull len=1;g==1&&len<=(1<<17);len*=2)
		{
			x=y;
			for(ull i=0;i<len;i++) y=Next(y);
			for(ull j=0;j<len&&g==1;j+=64)
			{
				ys=y;
				ull v=1;
				for(ull i=0;i<min<ull>(64,len-j);i++)
				{
					y=Next(y);
					v=Mul(v,x>y?x-y:y-x,n);
				}
				g=gcd(v,n);
			}
		}
		if(g==n)
		{
			do
			{
				ys=Next(ys);
				g=gcd(x>ys?x-ys:ys-x,n);
			}while(g==1);
		}
		if(g>1&&g<n) return g;
	}
	return 0;
}
ull Weight(ull r)
{
	return (S+r-2)/(r-1);
}
void Add(node &res,kind k,witness w)
{
	res.cnt[k]++;
	res.sum[k]+=Weight(w.r);
	if(!res.first[k].r||w.r<res.first[k].r) res.first[k]=w;
}
void Check(ull r,const vector<root> &ord,const vector<ull> &spec,node &res)
{
	ull n=r*r+r+1;
	if(r%3==1) n/=3;
	bool prime=IsPrime(n);
	ull a=n,b=0;
	if(prime) res.prime++;
	else
	{
		res.composite++;
		ull d=Factor(n);
		if(!d)
		{
			Add(res,unresolved,{r});
			return;
		}
		// All factors exceed D and n<D^3, so both factors are prime.
		a=min(d,n/d),b=max(d,n/d);
	}
	witness w;
	for(ull p:{a,b})
	{
		if(!p||p==w.p) continue;
		ull x=0,y=0;
		for(auto [d,c,e]:ord)
		{
			ull v=p%d;
			if(v==c||v==e)
			{
				x=d;
				break;
			}
		}
		if(!x) continue;
		if(p%4==1)
			for(ull d:spec)
				if(p%d==d-1)
				{
					y=d;
					break;
				}
		if(p%4==3||y)
		{
			witness ex{r,p,x,y};
			Add(res,p<r?back:(y?forwardSpec:forwardOrd),ex);
			if(prime) Add(res,forwardPrime,ex);
			return;
		}
		if(w.p)
		{
			Add(res,twoSpec,{r,w.p,w.a,0,p,x});
			return;
		}
		w={r,p,x};
	}
	if(w.p)
	{
		Add(res,cond,w);
		res.group[w.p]+=Weight(r);
	}
}
void Solve()
{
	vector<bool> vis(D+1);
	vector<ull> p,spec;
	vector<root> q,ord;
	for(ull i=3;i<=D;i+=2)
	{
		if(vis[i]) continue;
		if(i*i<=r) p.push_back(i);
		if(i*i<=D)
			for(ull j=i*i;j<=D;j+=2*i) vis[j]=1;
		if(i>3&&i<=lim) spec.push_back(i);
		if(i%3==1)
		{
			ull x=1;
			for(ull j=2;x==1;j++) x=Pows(j,(i-1)/3,i);
			root v{i,x,i-1-x};
			q.push_back(v);
			if(i<=lim) ord.push_back(v);
		}
	}
	atomic<ull> pos{(l+1)|1ULL};
	vector<node> ans(t);
	auto Work=[&](node &res)
	{
		vector<unsigned char> a(N);
		for(ull st=pos.fetch_add(2*N);st<=r;st=pos.fetch_add(2*N))
		{
			ull len=min<ull>(N,(r-st)/2+1),ed=st+2*(len-1);
			fill(a.begin(),a.begin()+len,0);
			for(ull x:p)
			{
				ull sq=x*x;
				if(sq>ed) break;
				ull f=max(sq,((st+x-1)/x)*x);
				if(!(f&1)) f+=x;
				for(ull j=(f-st)/2;j<len;j+=x) a[j]=1;
			}
			for(auto [x,y,z]:q)
			{
				ull rem=st%x;
				for(ull v:{y,z})
				{
					ull f=(v+x-rem)%x;
					if(f&1) f+=x;
					for(ull j=f/2;j<len;j+=x) a[j]=1;
				}
			}
			for(ull i=0;i<len;i++)
			{
				if(a[i]) continue;
				ull x=st+2*i;
				Add(res,base,{x});
				Check(x,ord,spec,res);
			}
		}
	};
	vector<thread> th;
	for(ull i=1;i<t;i++) th.emplace_back(Work,ref(ans[i]));
	Work(ans[0]);
	for(auto &x:th) x.join();
	node res;
	for(const auto &x:ans)
	{
		res.prime+=x.prime;
		res.composite+=x.composite;
		for(auto [p,s]:x.group) res.group[p]+=s;
		for(int i=0;i<K;i++)
		{
			res.cnt[i]+=x.cnt[i];
			res.sum[i]+=x.sum[i];
			if(x.first[i].r&&(!res.first[i].r||x.first[i].r<res.first[i].r))
				res.first[i]=x.first[i];
		}
	}
	ull pmax=0,smax=0;
	for(auto [p,s]:res.group)
		if(s>smax||(s==smax&&p<pmax)) pmax=p,smax=s;
	cout<<"special_groups\tmax_group_scaled_sum\tmax_group_prime\n";
	cout<<res.group.size()<<'\t'<<smax<<'\t'<<pmax<<'\n';
	cout<<"lower\tupper\tsecond_limit\tprime_quotients\tcomposite_quotients\n";
	cout<<l<<'\t'<<r<<'\t'<<lim<<'\t'<<res.prime<<'\t'<<res.composite<<'\n';
	cout<<"category\tcount\tscaled_sum\tr\tell\texcluded\tspecial\tother\tother_excluded\n";
	for(int i=0;i<K;i++)
	{
		auto [r,p,a,b,q,c]=res.first[i];
		cout<<name[i]<<'\t'<<res.cnt[i]<<'\t'<<res.sum[i]<<'\t'<<r<<'\t'
			<<p<<'\t'<<a<<'\t'<<b<<'\t'<<q<<'\t'<<c<<'\n';
	}
}
signed main(signed argc,char** argv)
{
	ios::sync_with_stdio(false);
	cin.tie(0);cout.tie(0);
	ull* a[]={&r,&t,&l,&lim};
	if(argc>5)
	{
		cerr<<"Usage: chain [upper [threads [lower [second_limit]]]]\n";
		return 1;
	}
	for(int i=1;i<argc;i++)
	{
		string_view s(argv[i]);
		auto x=from_chars(s.data(),s.data()+s.size(),*a[i-1]);
		if(x.ec!=errc{}||x.ptr!=s.data()+s.size())
		{
			cerr<<"Invalid argument.\n";
			return 1;
		}
	}
	if(l<D||r<=l||r>MAXN||t<1||t>128||lim<7||lim>D)
	{
		cerr<<"Invalid range or sieve parameters.\n";
		return 1;
	}
	Solve();
	return 0;
}
