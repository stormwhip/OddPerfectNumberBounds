# Odd Perfect Number Bounds — 受限奇完全数指数下界

沿用 Ochem–Zelinsky 的结构归约与已排除的小素数，通过分段同余筛、97 项容斥多项式以及显式素数计数误差界，对例外偶指数给出以下数值证书：

```math
q > 220000000000000000000 = 2.20 \times 10^{20}
```

## 项目文件

- [prefix.cpp](prefix.cpp)：C++17 多线程分段同余筛，计算候选素数数量及倒数和上界。
- [sieve.tsv](sieve.tsv)：4943 个连续筛选区间的输出数据。
- [verify.py](verify.py)：独立 Python 证明核验入口，检查 97 项容斥多项式、数据表汇总、显式误差系数及最终有理数不等式。
- [ResearchReport.md](ResearchReport.md)：数学推导、实验参数和结果。
- [DeterministicResearchPlan.md](DeterministicResearchPlan.md)：确定性数学与算法研究计划。

## 复现

```bash
python verify.py
g++ -std=c++17 -O3 -pthread prefix.cpp -o prefix
./prefix 50290000000000 16 50280000000000 3000000
```

`verify.py` 使用 Python 标准库，无需第三方 API；预期最后输出 `CERTIFIED: q > 220000000000000000000`。C++ 的演示区间应输出 `36314737` 与 `722196535949`。

验证器核验已保存的全部 4943 行筛选记录及其汇总；完整有限筛独立重算和解析误差定理的人工复核仍是投稿前的必要工作。

基础结构结果与因子树排除来自 Pascal Ochem、Joshua Zelinsky：[*On odd perfect numbers with exactly one even exponent greater than 2*](https://arxiv.org/abs/2607.19746)。
