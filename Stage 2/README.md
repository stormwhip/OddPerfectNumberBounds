# Stage 2 — 学习引导的容斥多项式搜索

本阶段研究机器学习能否在**可严格验证的筛函数搜索问题**中优于确定性算法。已完成的 LightGBM 实验属于探索性机器学习研究；扩展 LP 与约束生成实验属于本阶段的**非 ML 对照基线**。本阶段暂无新的可认证指数下界，当前最强严格结果见 [Stage 1](../Stage%201/)。

## 核心材料

- [LearningGuidedSievePlan.md](LearningGuidedSievePlan.md)：问题定义、研究方向、方法与实验计划，含前三轮探索记录。
- [PilotReport.md](experiments/PilotReport.md)：LightGBM 与确定性规则、随机策略比较，模型未超过确定性规则。
- [ExtendedSearchReport.md](experiments/ExtendedSearchReport.md)：12、14、15 素数、最高四次交集的完整 LP 实验。
- [ConstraintGenerationReport.md](experiments/ConstraintGenerationReport.md)：约束生成基线，在 16 素数问题上用 339 条活动约束代替 65535 条完整约束的求解，并核验所有布尔状态。

## 代码

```bash
# 在仓库根目录执行；Stage 2 需要 numpy/scipy，ML 实验另需 scikit-learn/lightgbm
python "Stage 2/experiments/lp_full_search.py" 12 100000
python "Stage 2/experiments/learning_guided_pilot.py" --seed 42
python "Stage 2/experiments/extended_lp_search.py" --primes 15 --degree 4 --maxmod 100000
python "Stage 2/experiments/constraint_generation.py" --primes 16
```

研究数据：[extended_search_summary.json](experiments/extended_search_summary.json)、[constraint_generation_summary.json](experiments/constraint_generation_summary.json)。以上 LP 数值实验均不能取代 Stage 1 的严格有理数证书。

## 下一步

现有低模数、低阶单项式池已被确定性算法有效处理。下一步优先构造**需要多轮列/约束生成的困难搜索问题**，再检验 LightGBM 或其他轻量模型能否优于纯 LP、贪心与确定性列生成。全部实验只使用本地开源模型与工具。
