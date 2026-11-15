# 受限奇完全数的指数下界

2026 年 11 月 16 日更新

## 摘要

设 $N$ 为奇完全数，且恰有一个偶指数大于 $2$。记这个素数幂为 $n^q$。沿用 Ochem 与 Zelinsky 对 $n=3$ 的归约、对 $3\,000\,000$ 以下素数的排除，以及引理 8，对小素数执行分段同余筛，对大素数用容斥法与等差数列素数计数的显式误差界估计。筛到 $B=50\,290\,000\,000\,000$ 后，得到

$$
\boxed{q>7.705\times10^{19}}.
$$

这个下界约为原论文 $2.3\times10^{10}$ 的 $3.35\times10^{9}$ 倍。有限筛和最后的不等式分别用整数和有理数计算。原论文的因子树和引理是本报告的输入。[原论文](https://arxiv.org/abs/2607.19746) · [原代码及数据](https://www.lirmm.fr/~ochem/opn/restrained/)

## 1. 把完全数条件化为倒数和

按照原论文第 4 节，写成 $N=3^q p^e\prod_{r\in R}r^2$，其中 $p^e$ 是 Euler 特殊分量，$R$ 是指数为 $2$ 的普通素因子集合。引理 8 给出 $q\ge \omega(N)/2-2$，因此 $|R|\le 2q+2$。已经排除的普通素数和特殊素数都可用于约束后面的筛选。[原论文第 4 节](https://arxiv.org/abs/2607.19746)

由于 $\sigma(3^q)/3^q<3/2$，以及任意素数幂 $u^a$ 满足

$$
\log\frac{\sigma(u^a)}{u^a}
<\log\left(1+\frac1{u-1}\right)<\frac1{u-1},
$$

由 $\sigma(N)=2N$ 得到必要条件

$$
\log(4/3)<\frac1{2\,999\,999}+\sum_{r\in R}\frac1{r-1}.
\tag{1}
$$

这里用 $p>3\,000\,000$ 估计特殊分量。选定 $Z>B$，其中 $B$ 是有限筛终点。给每个普通素因子预留 $1/(Z-1)$，再用

$$
w_Z(r)=\frac1{r-1}-\frac1{Z-1}\qquad(r\le Z)
$$

补足较小素数的贡献。于是只需对候选素数的 $w_Z$ 求和，再加上 $(2q+2)/(Z-1)$。这样可以计入每个小候选的实际大小。

## 2. 分段同余筛

若 $r\in R$，则 $r^2+r+1=\sigma(r^2)$ 的每个素因子都整除 $N$。对已排除的素数 $\ell\equiv1\pmod3$，$x^2+x+1\bmod\ell$ 的两个根就是模 $\ell$ 的非平凡三次单位根。程序先用埃拉托斯特尼筛去合数，再标记这两个剩余类。对每个剩余的素数累加 $\lceil 10^{18}/(r-1)\rceil$，向上取整保证总和是上界。[筛选程序](QBoundPrefix.cpp)

同一行的区间使用表内指定上限内的已排除素数。各段均包含其范围内所有可能的普通素因子。筛选程序使用 128 位整数累加，转换为 64 位无符号整数输出前检查数值范围。逐段记录及其筛除参数见数据表。[分段数据](QBoundSieveRuns.tsv)

| 完整区间 | 筛除所用素数的上限 | 候选数 | 倒数和上界，乘 $10^{18}$ |
|---|---:|---:|---:|
| $(3\,000\,000,200\,000\,000\,000]$ | 3 000 000 | 917 142 659 | 62 696 136 715 648 189 |
| $(200\,000\,000\,000,600\,000\,000\,000]$ | 100 000 | 2 225 037 956 | 6 133 899 512 944 946 |
| $(600\,000\,000\,000,50\,290\,000\,000\,000]$ | 3 000 000 | 186 184 186 178 | 17 324 085 008 638 533 |
| 合计 |  | 189 326 366 793 | 86 154 121 237 231 668 |

原 `Loopless.cpp` 在 64 位无符号整数里直接计算 `r * r + r + 1`，超过 $2^{32}$ 后会发生回绕。分段筛只计算模 $\ell$ 的根与各段起点，避开了这处溢出。原始数据 `Tree.txt` 的根覆盖了从 $11$ 到 $3\,000\,000$ 的全部 $216\,812$ 个素数，记录中列出的整除关系也已核对。小素数排除本身沿用论文及其因子树证明。[原筛选程序](OriginalWork/Source/Loopless.cpp) · [原因子树](OriginalWork/Source/Tree.txt)

## 3. 尾部界

对任意素数 $r>B$，令 $I_\ell(r)$ 表示 $\ell\mid r^2+r+1$。取

$$
H=\{13,19,31,37,43,61,67\},\qquad k(r)=\sum_{\ell\in H}I_\ell(r).
$$

通过这八个小素数筛选的指示函数至多为

$$
G(r)=(1-I_7(r))\left(1-k(r)+\binom{k(r)}2\right).
\tag{2}
$$

当 $I_7(r)=1$ 时，两边都为零。当 $I_7(r)=0$ 且 $k=0$ 时，两边都为一。$k\ge1$ 时，括号内为 $(k-1)(k-2)/2$，大于或等于零。这给出了逐项上界。

展开式 (2)，每一项最多涉及三个不同素数。若模数 $m$ 是这些素数的乘积，对应交集由 $2^j$ 个既约剩余类构成，其中 $j$ 是素数的个数。最大模数为 $7\cdot61\cdot67=28\,609$。对非空交集使用 Bennett 等人的定理 1.3，对空交集使用式 (1.15)。合并主项并对各误差取绝对值，可得

$$
\left|\sum_{\substack{r\le t\\r\ {\rm prime}}}G(r)-\delta\operatorname{Li}(t)\right|
\le C\frac{t}{\log^2 t},
\qquad
\delta=\frac{370141}{935550},\quad C=\frac{569}{840}.
\tag{3}
$$

式 (3) 中的模数范围保证 $t\ge8\times10^9$ 时可使用这些显式估计。主项系数来自把每个 $I_\ell$ 替换为 $2/(\ell-1)$。误差系数按交集的剩余类数量计入，模数不超过 $10^4$ 时每类取 $1/840$，其余每类取 $1/160$。有理数脚本直接从集合 $H$ 计算这两个系数。[等差数列素数计数定理](https://arxiv.org/abs/1802.00085v3) · [数值核验脚本](QBoundVerify.py)

设 $P_G(t)=\sum_{r\le t,\ r\ {\rm prime}}G(r)$。从 $w_Z(r)=\int_r^Z dt/(t-1)^2$ 得到尾部权重上界

$$
\sum_{\substack{B<r\le Z\\r\in R}} w_Z(r)
\le \int_B^Z \frac{P_G(t)-P_G(B)}{(t-1)^2}\,dt.
$$

将式 (3) 代入，主项和误差分别用下式估计。记 $b=\log B,\ z=\log Z,\ \rho=B/(B-1)$。

$$
\begin{aligned}
M&=\delta\left[\rho\log(z/b)
-\frac{Z/z-B/b+(Z-B)/z^2}{Z-1}\right],\\
E&=C\left[\frac{\rho}{b^2}+\rho^2\left(\frac1b-\frac1z\right)\right].
\end{aligned}
$$

主项中的 $Z/z-B/b+(Z-B)/z^2$ 是对 $\int_B^Z dt/\log t$ 做一次分部积分得到的下界。误差项使用 $t/(t-1)\le\rho$。因此两项之和是尾部权重的上界。

## 4. 数值结论与复算

取 $B=50\,290\,000\,000\,000, Z=18\,909\,040\,000\,000\,000\,000\,000$。设 $S=0.086154121237231668$ 为有限筛的倒数和上界，$K_B=189\,326\,366\,793$ 为该段候选数。若 $q\le 77\,050\,000\,000\,000\,000\,000$，式 (1) 的右边至多为

$$
\frac1{2\,999\,999}+S+M+E+\frac{2q+2-K_B}{Z-1}
\le 0.287681254791.
$$

有理数核验给出 $\log(4/3)\ge 0.287682072451$，与式 (1) 矛盾。脚本对对数使用带余项界的级数，对最终两个十进制数分别向上和向下取整。[数值核验脚本](QBoundVerify.py)

用已保存的数据直接复核最后一步。

```powershell
python QBoundVerify.py 50290000000000 18909040000000000000000 77050000000000000000 189326366793 86154121237231668
```

从头复算时，编译仓库中的 [QBoundPrefix.cpp](QBoundPrefix.cpp)。该版本允许筛选终点至 $10^{14}$，覆盖本报告中的 $B=5.029\times10^{13}$。下面的命令会重新运行完整有限筛，逐行核对候选数及倒数和上界。

```powershell
g++ -std=c++17 -O3 -pthread QBoundPrefix.cpp -o QBoundPrefix.exe
if ($LASTEXITCODE -ne 0) { throw "Compilation failed" }

$rows = Import-Csv .\QBoundSieveRuns.tsv -Delimiter ([char]9)
$count = 0L
$scaledSum = 0L
foreach ($row in $rows) {
    $result = & .\QBoundPrefix.exe $row.upper_inclusive 16 $row.lower_exclusive $row.residue_limit
    if ($LASTEXITCODE -ne 0 -or $result.Count -ne 2 -or
        $result[0] -ne $row.survivors -or $result[1] -ne $row.scaled_reciprocal_upper) {
        throw "Sieve result mismatch for endpoint $($row.upper_inclusive)"
    }
    $count += [long]$result[0]
    $scaledSum += [long]$result[1]
}
python QBoundVerify.py 50290000000000 18909040000000000000000 77050000000000000000 $count $scaledSum
```

## 5. 后续研究线索

进一步提升尾部界的方向，是调整式 (2) 的小素数组合，并同时比较主项系数 $\delta$ 和总误差系数 $C$。模数增大会增加每个剩余类的误差上界，因此两个系数都参与比较。

还可以用排除素数约束 $q+1$ 的因子。若 $\ell>3$ 已被排除，则 $\operatorname{ord}_\ell(3)\mid q+1$ 会迫使 $\ell\mid\sigma(3^q)$。例如 $\operatorname{ord}_{11}(3)=5$，所以 $5\nmid q+1$。这条线索可用于在数值下界之上排除特定的指数。

## 附录　问题 9 的一个解

赵睿泽用 Pell 方程构造得到素数 $x=47116128896261596331524418282573$。直接计算可得

$$
x^2+x+1=7^3\cdot43^2\cdot3416683^2\cdot12353161^2\cdot1401752363863729^2.
$$

右侧各因子均为素数，指数均至少为 $2$，所以该例回答了原论文问题 9。由于 $7$ 已在 $n=3$ 的情形中排除，这个 $x$ 会被同余筛排除。
