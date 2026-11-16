# Odd Perfect Number Bounds — 受限奇完全数指数下界

**本阶段已得到的结果：** 对受限奇完全数（restrained odd perfect numbers），在沿用 Ochem–Zelinsky 的结构归约及小素数排除条件下，得到严格数值证明

```math
q > 220000000000000000000 = 2.20 × 10^20
```

本仓库包含受限奇完全数指数下界研究的数学论证、有限筛程序、数据、精确证书，以及英文 LaTeX 论文初稿。

## 论文初稿

- [paper/main.tex](paper/main.tex)：仅围绕指数下界主定理的英文 LaTeX 初稿。
- [paper/AuditChecklist.md](paper/AuditChecklist.md)：证明严谨性与计算复现核验清单。

## 文件与依赖

- [ResearchReport.md](ResearchReport.md)：正文数学推导，97 项容斥多项式、显式素数计数尾部估计、有限筛。
- [QBoundVerifyAPRefined.py](QBoundVerifyAPRefined.py)：当前 2.20×10^20 下界的精确验证。
- [DeterministicResearchPlan.md](DeterministicResearchPlan.md)：确定性算法改进研究路线。
- [QBoundVerifyOptimized.py](QBoundVerifyOptimized.py)：前一版 2.193×10^20 下界的**精确有理数验证**，检查 4096 个布尔状态、97 项、4943 段汇总。
- [QBoundVerify.py](QBoundVerify.py)：最初的 7.705×10^19 下界及供优化验证器调用的有理数计算。
- [QBoundPrefix.cpp](QBoundPrefix.cpp)：C++17 分段筛程序。
- [QBoundSieveRuns.tsv](QBoundSieveRuns.tsv)：保存的 4943 段筛选数据。
- [OriginalWork](OriginalWork/)：原作者提供的论文及伴随程序、因子树数据，保留原文件内容，并在总 README 中署名。

## 复现

从仓库根目录执行：

```bash
python QBoundVerifyAPRefined.py
python QBoundVerifyOptimized.py
python QBoundVerify.py 50290000000000 18909040000000000000000 77050000000000000000 189326366793 86154121237231668
g++ -std=c++17 -O3 -pthread QBoundPrefix.cpp -o QBoundPrefix
./QBoundPrefix 50290000000000 16 50280000000000 3000000
```

首个命令应给出 `CERTIFIED: q > 220000000000000000000`，第二个命令验证前一版下界。C++ 最后一个演示区间应输出 `36314737` 与 `722196535949`。完整有限筛独立重跑仍属于正式投稿前的必要工作。

论文与程序来源：Pascal Ochem、Joshua Zelinsky，[*On odd perfect numbers with exactly one even exponent greater than 2*](https://arxiv.org/abs/2607.19746)。
