# 首轮实验：学习引导的筛函数搜索（2026-11-16）

## 结果摘要

**本轮已经实际训练并测试 LightGBM；没有产生新的可认证数论下界。**

当前记录的严格结果仍为 q > 2.193 × 10^20，证书由仓库的 QBoundVerifyOptimized.py 验证。

### A. 纯线性规划基线

在固定 B=50,290,000,000,000、Z=57,433,500,000,000,000,000,000 下，对候选多项式系数建立浮点线性规划：最小化 A·delta + D·C，且对每个非零布尔向量 z 都有 G(z)≥0，G(0)=1；系数绝对值以拆分变量处理。数值实验采用 SciPy/HiGHS，成功返回最优状态。

| 素数集合 | 模数上限 | 合法单项式 | 求解时间（本地单次） | 结果 |
| --- | ---: | ---: | ---: | --- |
| 前 12 个素数 | 10,000 | 99 | 0.55 秒 | 当前 97 项证书 |
| 前 12 个素数 | 100,000 | 211 | 1.09 秒 | 当前 97 项证书 |
| 前 13 个素数（增加 109） | 100,000 | 244 | 4.20 秒 | 当前 97 项证书 |

三个实验均得到 delta=0.3829121454979...、C=0.515476190476...，实际非零项为原证书的 97 项。因固定 Z 和允许项集合，这说明在**上述特定线性规划模型**下已达到最优；不代表一般筛法或全局 q 下界达到极限。

### B. LightGBM 与确定性/随机排序的受控对照

可复现脚本：[learning_guided_pilot.py](learning_guided_pilot.py)。

- 数据：每次从 12 个已排除素数中选择含 7 的 8–11 元子集，固定随机种子；每次用完整 LP 产生训练标签（哪些三次项在最优解中非零）。
- 划分：每个种子 60 个训练实例、20 个**未见过的素数子集**测试实例；没有让同一子集的不同状态跨越训练/测试分区。
- 模型：本地 LightGBM 分类器，70 棵树；预测三次单项式是否属于最优 LP 解。
- 比较：在每个测试实例中只允许选择 5 个三次候选项，其余常数、一次、二次单项式全部保留；再用同一完整布尔约束的 LP 求解。与包含 7 的规则+净收益排序、随机选项和全部合法列 LP 作比较。
- 指标：相对“未添加三次项”的尾部目标减小量。**该指标是浮点优化目标，不能解释为已经证明的正 q 下界。**

每个随机种子的得分为“取得完整 LP 所能实现的目标改善比例”：

| 随机种子 | LightGBM | 确定性规则 | 随机选项 |
| --- | ---: | ---: | ---: |
| 42 | 71.89% | 75.04% | 55.11% |
| 123 | 71.83% | 74.00% | 53.95% |
| 2026 | 73.18% | 74.75% | 52.73% |
| 2027 | 73.90% | 76.45% | 56.18% |
| **四组平均** | **72.70%** | **75.06%** | **54.49%** |

每组测试中，LightGBM 对最优 LP 选中项的分类 ROC-AUC 和 Average Precision 都为 1.0。这一分类任务**非常简单**：最优解的三次项共享素数 7，而这一特征已经显式提供给模型。分类准确度不足以体现优于数学规则的搜索能力。四组中，模型在每组 20 个测试实例里均未在目标改善量上超过确定性规则。

### C. 研究判断

当前小规模的合法单项式空间已经可以由直接 LP 在数秒内求解。机器学习在该实验中优于随机排序，但逊于廉价的手工排序规则；直接将 ML 作为本文主贡献的论据不足。

下一轮优先转向**大规模、列生成昂贵的实例**，例如更多筛选素数、更高阶合法交集、模数成本分段、联合优化 B/Z 与学习策略，并与 reduced-cost 列生成进行严格的耗时及最优性对照。单项式允许的模数必须逐一满足 Bennett–Martin–O'Bryant–Rechnitzer 的显式定理适用条件。所有浮点候选都须由独立有理数证书验证。

## 复现

安装 Python 3.11+ 的 numpy、scipy、scikit-learn、lightgbm；在仓库根目录依次执行：

```bash
python research/experiments/learning_guided_pilot.py --train 60 --test 20 --budget 5 --seed 42 --output research/experiments/seed42.json
python research/experiments/learning_guided_pilot.py --train 60 --test 20 --budget 5 --seed 123 --output research/experiments/seed123.json
python research/experiments/learning_guided_pilot.py --train 60 --test 20 --budget 5 --seed 2026 --output research/experiments/seed2026.json
python research/experiments/learning_guided_pilot.py --train 60 --test 20 --budget 5 --seed 2027 --output research/experiments/seed2027.json
python QBoundVerifyOptimized.py
```

环境：本地 CPU、SciPy/HiGHS 与 LightGBM；未使用网络推理、第三方 API 密钥或付费模型。

**论文状态**：继续沿用原已验证的 q 下界；机器学习部分属于探索性负结果，不作为新的数学定理。
