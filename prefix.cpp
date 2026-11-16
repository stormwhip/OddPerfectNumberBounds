/*
Build: g++ -std=c++17 -O2 -pthread prefix.cpp -o prefix
Usage: ./prefix [upper [threads [lower [residue_limit]]]]
Output: surviving candidate count and reciprocal-sum upper bound multiplied by 10^18.
*/
#include<bits/stdc++.h>
#define int long long
using namespace std;
using ull=unsigned long long;
using u128=__uint128_t;
const int N=1<<19;
const ull S=1000000000000000000ULL,MAXN=100000000000000ULL;
struct node
{
	ull cnt=0;
	u128 sum=0;
};
ull l=3000000,r=10000000000ULL,t=1,lim=100000;
int Pows(int a,int b,int mod)
{
	int ans=1;
	while(b)
	{
		if(b&1) ans=ans*a%mod;
		a=a*a%mod;
		b>>=1;
	}
	return ans;
}
void Solve()
{
	int n=lim;
	while((ull)n*n<r) n++;
	vector<bool> vis(n+1);
	vector<int> p;
	vector<pair<int,int>> q;
	for(int i=3;i<=n;i+=2)
	{
		if(vis[i]) continue;
		if((ull)i*i<=r) p.push_back(i);
		if(i<=n/i)
			for(int j=i*i;j<=n;j+=2*i) vis[j]=1;
		if((ull)i<=lim&&i%3==1)
		{
			int x=1;
			for(int j=2;x==1;j++) x=Pows(j,(i-1)/3,i);
			q.push_back({i,x});
			q.push_back({i,i-1-x});
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
			for(int x:p)
			{
				ull sq=(ull)x*x;
				if(sq>ed) break;
				ull f=max(sq,((st+x-1)/x)*x);
				if(!(f&1)) f+=x;
				for(ull j=(f-st)/2;j<len;j+=x) a[j]=1;
			}
			for(auto [x,y]:q)
			{
				ull f=st+(y+x-st%x)%x;
				if(!(f&1)) f+=x;
				for(ull j=(f-st)/2;j<len;j+=x) a[j]=1;
			}
			for(ull i=0;i<len;i++)
			{
				if(a[i]) continue;
				res.cnt++;
				ull d=st+2*i-1;
				res.sum+=(S+d-1)/d;
			}
		}
	};
	vector<thread> th;
	for(int i=1;i<t;i++) th.emplace_back(Work,ref(ans[i]));
	Work(ans[0]);
	for(auto &x:th) x.join();
	node res;
	for(auto x:ans) res.cnt+=x.cnt,res.sum+=x.sum;
	if(res.sum>numeric_limits<ull>::max())
	{
		cerr<<"Scaled reciprocal sum exceeds the output range.\n";
		exit(1);
	}
	cout<<res.cnt<<"\n"<<(ull)res.sum<<"\n";
}
signed main(signed argc,char** argv)
{
	ios::sync_with_stdio(false);
	cin.tie(0);cout.tie(0);
	ull* a[]={&r,&t,&l,&lim};
	if(argc>5)
	{
		cerr<<"Usage: prefix [upper [threads [lower [residue_limit]]]]\n";
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
	if(l<3000000||r<=l||r>MAXN||t<1||t>128||lim>3000000)
	{
		cerr<<"Invalid range or sieve parameters.\n";
		return 1;
	}
	Solve();
	return 0;
}
