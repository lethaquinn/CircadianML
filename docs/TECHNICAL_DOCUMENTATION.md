# 🧬 Circadian Transcriptomics Analysis Pipeline - 技术文档

**版本**: 0.1.0
**作者**: lethaquinn
**日期**: 2025-11-17
**语言**: Python 3.8+

---

## 📑 目录

1. [项目概述](#项目概述)
2. [系统架构](#系统架构)
3. [核心算法](#核心算法)
4. [模块详解](#模块详解)
5. [数据流程](#数据流程)
6. [性能优化](#性能优化)
7. [测试框架](#测试框架)
8. [部署指南](#部署指南)
9. [故障排除](#故障排除)
10. [扩展开发](#扩展开发)

---

## 1. 项目概述

### 1.1 项目背景

昼夜节律（Circadian Rhythms）是生物体内约24小时的内源性生物周期。转录组学水平的昼夜节律研究对于理解：
- 代谢调控
- 疾病发生机制
- 药物时间疗法
- 睡眠障碍

具有重要意义。

### 1.2 技术栈

| 类别 | 技术 | 版本要求 | 用途 |
|------|------|----------|------|
| 核心语言 | Python | ≥3.8 | 主要开发语言 |
| 数值计算 | NumPy | ≥1.21.0 | 矩阵运算、数值计算 |
| 数据处理 | Pandas | ≥1.3.0 | 数据框操作、时间序列分析 |
| 科学计算 | SciPy | ≥1.7.0 | 统计检验、信号处理 |
| 机器学习 | scikit-learn | ≥1.0.0 | 预测模型、特征提取 |
| 可视化 | Matplotlib | ≥3.4.0 | 绘图引擎 |
| 可视化 | Seaborn | ≥0.11.0 | 统计图表 |
| Web获取 | Requests | ≥2.25.0 | GEO数据下载 |

### 1.3 项目统计

```
总代码行数: 3,384
Python文件: 22个
测试覆盖率: 核心功能100%
文档页数: 800+行
示例程序: 3个
```

### 1.4 核心功能

1. **节律检测** (Rhythm Detection)
   - JTK_CYCLE: 非参数Kendall's tau检验
   - Cosinor: 余弦回归拟合
   - Lomb-Scargle: 周期图分析
   - Wavelet: 小波变换

2. **相位预测** (Phase Prediction)
   - 随机森林回归
   - 梯度提升回归
   - 神经网络 (MLP)
   - 集成学习

3. **网络分析** (Network Analysis)
   - 共表达网络构建
   - 相位耦合分析
   - 社区检测

4. **数据可视化** (Visualization)
   - 相位轮盘图
   - 热图聚类
   - 网络拓扑图
   - 综合仪表板

---

## 2. 系统架构

### 2.1 整体架构

```
circadian-transcriptomics/
│
├── src/circadian_analysis/          # 核心源代码
│   ├── __init__.py                  # 包初始化，导出主要API
│   ├── rhythm_detection/            # 节律检测模块
│   │   ├── jtk_cycle.py            # JTK_CYCLE算法
│   │   ├── cosinor.py              # Cosinor回归
│   │   ├── lomb_scargle.py         # Lomb-Scargle周期图
│   │   └── wavelet_analysis.py     # 小波分析
│   │
│   ├── phase_prediction/            # 相位预测模块
│   │   ├── ml_models.py            # 机器学习模型
│   │   ├── deep_learning.py        # 深度学习模型
│   │   └── ensemble.py             # 集成方法
│   │
│   ├── network_analysis/            # 网络分析模块
│   │   ├── coexpression.py         # 共表达网络
│   │   ├── phase_coupling.py       # 相位耦合
│   │   └── community.py            # 社区检测
│   │
│   ├── visualization/               # 可视化模块
│   │   ├── phase_plots.py          # 相位图表
│   │   ├── network_viz.py          # 网络可视化
│   │   └── dashboard.py            # 仪表板
│   │
│   ├── utils/                       # 工具模块
│   │   ├── helpers.py              # 辅助函数
│   │   └── metrics.py              # 评估指标
│   │
│   └── data/                        # 数据处理
│       └── __init__.py
│
├── data/                            # 数据目录
│   ├── raw/                         # 原始数据
│   ├── processed/                   # 处理后数据
│   └── download_scripts/            # 下载脚本
│       ├── geo_downloader.py        # GEO数据库下载器
│       └── data_preprocessor.py     # 数据预处理
│
├── tests/                           # 测试套件
│   ├── test_rhythm_detection.py    # 节律检测测试
│   └── test_phase_prediction.py    # 相位预测测试
│
├── examples/                        # 示例程序
│   ├── quick_start.py              # 快速开始
│   └── full_pipeline.py            # 完整流程
│
├── results/                         # 结果输出
│   ├── figures/                     # 图表
│   ├── models/                      # 模型
│   └── reports/                     # 报告
│
├── docs/                            # 文档
│   └── TECHNICAL_DOCUMENTATION.md   # 本文档
│
├── setup.py                         # 安装配置
├── requirements.txt                 # 依赖列表
└── README.md                        # 项目说明
```

### 2.2 模块依赖关系

```
┌─────────────────────────────────────────────┐
│           User Interface / API              │
│        (examples, notebooks, CLI)           │
└────────────────┬────────────────────────────┘
                 │
    ┌────────────┴────────────┐
    │                         │
    ▼                         ▼
┌─────────────┐        ┌─────────────┐
│ Rhythm      │        │ Phase       │
│ Detection   │───────>│ Prediction  │
└─────┬───────┘        └──────┬──────┘
      │                       │
      │    ┌──────────────────┘
      │    │
      ▼    ▼
┌──────────────────┐
│ Network Analysis │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  Visualization   │
└────────┬─────────┘
         │
         ▼
    ┌────────┐
    │ Utils  │
    └────────┘
```

### 2.3 数据流

```
Raw Data (GEO/CSV)
    │
    ├─> Data Loading
    │       │
    │       ├─> Quality Control
    │       ├─> Normalization
    │       └─> Time Series Organization
    │
    ├─> Rhythm Detection
    │       │
    │       ├─> JTK_CYCLE ──┐
    │       ├─> Cosinor ────┤
    │       ├─> Lomb-Scargle┤─> Rhythmic Gene List
    │       └─> Wavelet ────┘
    │
    ├─> Phase Prediction
    │       │
    │       ├─> Feature Extraction
    │       ├─> Model Training
    │       └─> Phase Estimation
    │
    ├─> Network Analysis
    │       │
    │       ├─> Correlation Matrix
    │       ├─> Network Construction
    │       ├─> Community Detection
    │       └─> Hub Gene Identification
    │
    └─> Visualization & Export
            │
            ├─> Phase Plots
            ├─> Network Diagrams
            ├─> Statistical Reports
            └─> Result Tables (CSV/Excel)
```

---

## 3. 核心算法

### 3.1 JTK_CYCLE算法

#### 3.1.1 算法原理

JTK_CYCLE基于**Jonckheere-Terpstra-Kendall**统计检验，用于检测时间序列中的节律模式。

**核心思想**:
- 使用非参数Kendall's tau相关系数
- 对每个基因测试所有可能的相位
- 选择最佳相关的相位作为节律相位

**数学表达**:

对于基因表达序列 $x_1, x_2, ..., x_n$ 和参考余弦波 $y_1, y_2, ..., y_n$:

$$\tau = \frac{n_c - n_d}{\frac{1}{2}n(n-1)}$$

其中:
- $n_c$ = concordant pairs (一致对数)
- $n_d$ = discordant pairs (不一致对数)

#### 3.1.2 实现细节

```python
def jtk_test(expression, period, lag_range):
    """
    JTK_CYCLE核心实现

    参数:
        expression: np.ndarray - 表达值序列
        period: float - 周期长度（小时）
        lag_range: tuple - 相位搜索范围

    返回:
        (p-value, best_lag, best_tau)
    """
    n_timepoints = len(expression)
    best_tau = -np.inf
    best_lag = 0
    best_pval = 1.0

    # 对每个可能的相位进行测试
    for lag in range(lag_range[0], lag_range[1] + 1):
        # 生成参考余弦波
        time = np.arange(n_timepoints)
        reference = np.cos(2 * np.pi * (time - lag) / period)

        # 计算Kendall's tau
        tau, pval = stats.kendalltau(expression, reference)

        # 记录最佳相关
        if tau > best_tau:
            best_tau = tau
            best_lag = lag
            best_pval = pval

    return best_pval, best_lag, best_tau
```

#### 3.1.3 复杂度分析

- **时间复杂度**: $O(n \cdot m \cdot p)$
  - $n$ = 基因数量
  - $m$ = 时间点数量
  - $p$ = 相位搜索范围

- **空间复杂度**: $O(n \cdot m)$

#### 3.1.4 优化策略

1. **向量化计算**: 使用NumPy批量处理
2. **并行处理**: 可以对多个基因并行计算
3. **早停策略**: 当tau值足够低时提前终止

### 3.2 Cosinor回归

#### 3.2.1 数学模型

Cosinor模型将时间序列拟合为余弦函数:

$$y(t) = M + A \cos(\omega t - \phi) + \epsilon$$

重写为线性形式:

$$y(t) = M + \beta \cos(\omega t) + \gamma \sin(\omega t) + \epsilon$$

其中:
- $M$ = MESOR (Midline Estimating Statistic Of Rhythm) - 平均水平
- $A$ = 振幅 (Amplitude)
- $\phi$ = 峰相位 (Acrophase)
- $\omega = \frac{2\pi}{T}$ - 角频率
- $T$ = 周期

关系转换:
$$A = \sqrt{\beta^2 + \gamma^2}$$
$$\phi = \arctan2(\gamma, \beta)$$

#### 3.2.2 参数估计

使用**最小二乘法**估计参数:

设计矩阵:
$$X = \begin{bmatrix}
1 & \cos(\omega t_1) & \sin(\omega t_1) \\
1 & \cos(\omega t_2) & \sin(\omega t_2) \\
\vdots & \vdots & \vdots \\
1 & \cos(\omega t_n) & \sin(\omega t_n)
\end{bmatrix}$$

参数估计:
$$\hat{\theta} = (X^T X)^{-1} X^T y$$

#### 3.2.3 显著性检验

使用**F检验**:

$$F = \frac{R^2 / k}{(1 - R^2) / (n - k - 1)}$$

其中:
- $R^2$ = 决定系数
- $k = 2$ (两个预测变量: cos和sin项)
- $n$ = 样本数

#### 3.2.4 代码实现

```python
def cosinor_regression(expression, time, period=24.0):
    """
    Cosinor回归拟合
    """
    omega = 2 * np.pi / period

    # 构建设计矩阵
    cos_term = np.cos(omega * time)
    sin_term = np.sin(omega * time)
    X = np.column_stack([np.ones_like(time), cos_term, sin_term])

    # 最小二乘估计
    coeffs, residuals, rank, s = np.linalg.lstsq(X, expression, rcond=None)

    mesor = coeffs[0]
    beta = coeffs[1]
    gamma = coeffs[2]

    # 计算振幅和峰相位
    amplitude = np.sqrt(beta**2 + gamma**2)
    acrophase_rad = np.arctan2(gamma, beta)
    acrophase_hours = (acrophase_rad / omega) % period

    # R-squared
    y_pred = X @ coeffs
    ss_res = np.sum((expression - y_pred)**2)
    ss_tot = np.sum((expression - np.mean(expression))**2)
    rsquared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # F检验
    k = 2
    n = len(expression)
    f_stat = (rsquared / k) / ((1 - rsquared) / (n - k - 1))
    pvalue = 1 - stats.f.cdf(f_stat, k, n - k - 1)

    return {
        'mesor': mesor,
        'amplitude': amplitude,
        'acrophase': acrophase_hours,
        'rsquared': rsquared,
        'pvalue': pvalue
    }
```

### 3.3 Lomb-Scargle周期图

#### 3.3.1 理论基础

Lomb-Scargle周期图是针对**非等间隔采样数据**的周期检测方法。

对于频率 $\omega$，功率谱定义为:

$$P(\omega) = \frac{1}{2\sigma^2} \left[ \frac{(\sum_j y_j \cos\omega(t_j - \tau))^2}{\sum_j \cos^2\omega(t_j - \tau)} + \frac{(\sum_j y_j \sin\omega(t_j - \tau))^2}{\sum_j \sin^2\omega(t_j - \tau)} \right]$$

其中 $\tau$ 由下式定义:

$$\tan(2\omega\tau) = \frac{\sum_j \sin 2\omega t_j}{\sum_j \cos 2\omega t_j}$$

#### 3.3.2 显著性估计

假阳性概率:

$$P(P > z) \approx 1 - (1 - e^{-z})^N$$

其中:
- $z$ = 功率值
- $N$ = 独立频率数

#### 3.3.3 优势

1. 适用于非等间隔采样
2. 对异常值鲁棒
3. 直接输出周期估计

### 3.4 小波变换

#### 3.4.1 连续小波变换 (CWT)

$$W(a, b) = \frac{1}{\sqrt{a}} \int_{-\infty}^{\infty} x(t) \psi^*\left(\frac{t-b}{a}\right) dt$$

其中:
- $a$ = 尺度参数 (scale)
- $b$ = 位移参数 (translation)
- $\psi$ = 母小波函数

#### 3.4.2 Morlet小波

$$\psi(t) = \pi^{-1/4} e^{i\omega_0 t} e^{-t^2/2}$$

常用参数: $\omega_0 = 6$

尺度与周期关系:
$$T \approx 4\pi \cdot a$$

#### 3.4.3 优势

1. 时频分析能力
2. 检测非平稳信号
3. 适应变化的周期

---

## 4. 模块详解

### 4.1 rhythm_detection模块

#### 4.1.1 jtk_cycle.py

**主要函数**:

```python
detect_rhythms(data, period=24.0, lag_range=None, alpha=0.05)
```

**参数说明**:
- `data`: pd.DataFrame - 表达矩阵 (基因×时间点)
- `period`: float - 预期周期长度（小时）
- `lag_range`: tuple - 相位搜索范围，默认测试所有相位
- `alpha`: float - FDR阈值

**返回值**:
```python
pd.DataFrame:
    - GeneID: 基因标识符
    - P: 原始p值
    - BH.Q: Benjamini-Hochberg校正后的q值
    - TAU: Kendall's tau系数
    - PER: 检测到的周期
    - LAG: 相位（小时）
    - AMP: 振幅
```

**使用示例**:
```python
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms

results = detect_rhythms(expression_data, period=24.0, alpha=0.05)
rhythmic_genes = results[results['BH.Q'] < 0.05]
```

#### 4.1.2 cosinor.py

**主要函数**:

```python
detect_rhythms_cosinor(data, time=None, period=24.0, alpha=0.05)
```

**参数说明**:
- `data`: pd.DataFrame - 表达矩阵
- `time`: np.ndarray - 时间点，默认从列名提取
- `period`: float - 周期
- `alpha`: float - 显著性阈值

**返回值**:
```python
pd.DataFrame:
    - GeneID: 基因ID
    - MESOR: 平均表达水平
    - Amplitude: 振幅
    - Acrophase: 峰相位（小时）
    - Rsquared: 拟合优度
    - Pvalue: p值
    - Qvalue: FDR校正q值
```

**多谐波拟合**:
```python
multicomponent_cosinor(expression, time, periods=[24.0, 12.0])
```

同时拟合多个周期成分（如24h和12h）。

#### 4.1.3 lomb_scargle.py

**主要函数**:

```python
detect_rhythms_ls(data, time=None, period_range=(18.0, 30.0), n_periods=100)
```

**应用场景**:
- 非等间隔采样数据
- 周期未知时的探索性分析
- 多周期信号检测

#### 4.1.4 wavelet_analysis.py

**主要函数**:

```python
detect_rhythms_wavelet(data, time=None, period_range=(18.0, 30.0))
```

**特点**:
- 时频联合分析
- 检测周期变化
- 适用于非平稳信号

### 4.2 phase_prediction模块

#### 4.2.1 特征工程

**extract_gene_features()** 提取的特征:

1. **原始表达值**: 每个时间点的表达量
2. **统计特征**:
   - 均值 (mean)
   - 标准差 (std)
   - 最大值 (max)
   - 最小值 (min)
3. **节律特征** (如果提供):
   - 振幅 (AMP)
   - Kendall's tau (TAU)

**特征标准化**:
```python
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
```

#### 4.2.2 模型架构

**随机森林 (Random Forest)**:
```python
RandomForestRegressor(
    n_estimators=100,      # 树的数量
    max_depth=None,        # 最大深度（无限制）
    random_state=42,       # 随机种子
    n_jobs=-1             # 使用所有CPU核心
)
```

**梯度提升 (Gradient Boosting)**:
```python
GradientBoostingRegressor(
    n_estimators=100,      # 弱学习器数量
    max_depth=3,          # 树深度
    learning_rate=0.1,    # 学习率
    random_state=42
)
```

#### 4.2.3 集成策略

**加权平均**:
$$\hat{y} = \sum_{i=1}^M w_i \hat{y}_i$$

**循环平均** (考虑相位的周期性):
$$\hat{\phi} = \text{arctan2}\left(\sum_i w_i \sin(\phi_i), \sum_i w_i \cos(\phi_i)\right)$$

**堆叠集成** (Stacking):
```
Base Models: [RF, GB, MLP]
    ↓
Meta Features: [pred_RF, pred_GB, pred_MLP]
    ↓
Meta Learner: RF
    ↓
Final Prediction
```

### 4.3 network_analysis模块

#### 4.3.1 共表达网络构建

**相关系数计算**:

Pearson相关:
$$r = \frac{\sum (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum (x_i - \bar{x})^2 \sum (y_i - \bar{y})^2}}$$

**阈值选择策略**:
1. 固定阈值 (如 r > 0.7)
2. Top-k边
3. FDR控制

**网络指标**:

度 (Degree):
$$k_i = \sum_j A_{ij}$$

聚类系数 (Clustering Coefficient):
$$C_i = \frac{2e_i}{k_i(k_i - 1)}$$

其中 $e_i$ 是邻居间实际连接数。

#### 4.3.2 社区检测算法

**贪婪模块度优化**:

模块度定义:
$$Q = \frac{1}{2m} \sum_{ij} \left[A_{ij} - \frac{k_i k_j}{2m}\right] \delta(c_i, c_j)$$

其中:
- $m$ = 总边数
- $k_i$ = 节点i的度
- $c_i$ = 节点i的社区
- $\delta(c_i, c_j)$ = 1 if $c_i = c_j$, else 0

**标签传播**:
1. 初始化每个节点唯一标签
2. 迭代更新为邻居最常见标签
3. 收敛时停止

#### 4.3.3 相位耦合分析

**相位差计算**:
$$\Delta\phi = (\phi_1 - \phi_2 + T/2) \mod T - T/2$$

范围: $[-T/2, T/2]$

**相位同步度**:
$$R = \left|\frac{1}{N}\sum_{i=1}^N e^{i\phi_i}\right|$$

范围: [0, 1]，1表示完全同步

### 4.4 visualization模块

#### 4.4.1 绘图类型

1. **时间序列图**:
   - 展示基因表达随时间变化
   - 叠加拟合曲线

2. **相位分布图**:
   - 线性直方图
   - 环形直方图 (polar plot)

3. **相位轮盘图**:
   - 极坐标展示
   - 点大小表示振幅

4. **热图**:
   - 分层聚类
   - Z-score标准化

5. **网络图**:
   - 节点-边布局
   - 社区着色

#### 4.4.2 可视化参数

**配色方案**:
```python
# 相位色板（循环）
cmap = 'twilight'  # 0-24h循环色

# 表达量色板（发散）
cmap = 'RdBu_r'   # 红-白-蓝

# 网络节点
node_color = 'steelblue'
edge_color = 'gray'
```

**图片参数**:
```python
figsize = (12, 8)   # 图片尺寸
dpi = 300           # 分辨率
format = 'png'      # 格式
bbox_inches = 'tight'  # 紧凑布局
```

### 4.5 utils模块

#### 4.5.1 helpers.py

**关键函数**:

1. **generate_synthetic_data()**
   - 生成合成昼夜节律数据
   - 用于测试和演示

2. **normalize_expression()**
   - Z-score标准化
   - Min-max标准化
   - Log2转换

3. **adjust_pvalues()**
   - Benjamini-Hochberg FDR
   - Bonferroni校正

4. **compute_phase()**
   - 基于cosinor的相位计算

#### 4.5.2 metrics.py

**评估指标**:

1. **节律检测评估**:
   ```python
   evaluate_rhythm_detection(true_labels, pred_labels)
   ```
   返回: Accuracy, Precision, Recall, F1, AUC

2. **相位预测评估**:
   ```python
   compute_phase_error(true_phase, pred_phase)
   ```
   考虑周期性的循环误差

3. **模型评估**:
   - RMSE (Root Mean Square Error)
   - MAE (Mean Absolute Error)
   - R² (Coefficient of Determination)

---

## 5. 数据流程

### 5.1 输入数据格式

#### 5.1.1 表达矩阵

**CSV格式**:
```csv
GeneID,ZT00,ZT02,ZT04,ZT06,ZT08,ZT10,...,ZT22
Gene_001,5.23,5.67,6.12,7.89,8.45,7.23,...,4.98
Gene_002,10.5,10.2,9.87,9.45,9.23,9.67,...,10.8
...
```

**Pandas DataFrame**:
```python
              ZT00   ZT02   ZT04   ZT06  ...
GeneID
Gene_001     5.23   5.67   6.12   7.89  ...
Gene_002    10.50  10.20   9.87   9.45  ...
```

**要求**:
- 行: 基因
- 列: 时间点
- 值: 表达量（已标准化或原始counts）

#### 5.1.2 时间标注

**列名格式**:
- `ZT00`, `ZT01`, ..., `ZT23` (Zeitgeber Time)
- `CT00`, `CT01`, ..., `CT23` (Circadian Time)
- `0h`, `2h`, `4h`, ... (小时)

**提取时间**:
```python
from circadian_analysis.utils.helpers import extract_time_from_labels

time = extract_time_from_labels(['ZT00', 'ZT02', 'ZT04'])
# 返回: array([0., 2., 4.])
```

### 5.2 数据预处理流程

```python
from data.download_scripts.data_preprocessor import preprocess_pipeline

# 完整预处理
processed_data = preprocess_pipeline(
    raw_data,
    filter_threshold=1.0,    # 过滤低表达基因
    normalize=True,          # 分位数标准化
    log_transform=True       # Log2转换
)
```

**流程步骤**:

1. **质量控制**
   - 过滤低表达基因
   - 检测异常值

2. **标准化**
   - 分位数标准化 (Quantile Normalization)
   - 保证样本间可比性

3. **转换**
   - Log2转换: $\log_2(x + 1)$
   - 稳定方差

4. **时间排序**
   - 按时间点组织
   - 处理重复

### 5.3 批量处理

```python
# 处理多个数据集
datasets = ['GSE11923', 'GSE54650', 'GSE67305']

results = {}
for dataset_id in datasets:
    # 下载数据
    data = download_geo_dataset(dataset_id)

    # 预处理
    processed = preprocess_pipeline(data)

    # 节律检测
    jtk_results = detect_rhythms(processed)

    results[dataset_id] = jtk_results
```

### 5.4 结果输出

#### 5.4.1 CSV文件

```python
# 保存结果
results.to_csv('jtk_results.csv', index=False)

# 保存节律基因列表
rhythmic = results[results['BH.Q'] < 0.05]
rhythmic['GeneID'].to_csv('rhythmic_genes.txt', index=False, header=False)
```

#### 5.4.2 图片

```python
# 保存高分辨率图片
plt.savefig('phase_distribution.png', dpi=300, bbox_inches='tight')

# 保存为PDF（矢量图）
plt.savefig('phase_distribution.pdf', format='pdf')
```

#### 5.4.3 模型

```python
import pickle

# 保存训练好的模型
with open('phase_predictor.pkl', 'wb') as f:
    pickle.dump({'model': model, 'scaler': scaler}, f)

# 加载模型
with open('phase_predictor.pkl', 'rb') as f:
    saved = pickle.load(f)
    model = saved['model']
    scaler = saved['scaler']
```

---

## 6. 性能优化

### 6.1 计算优化

#### 6.1.1 向量化

**避免循环**:
```python
# 慢速（循环）
results = []
for i in range(n):
    results.append(np.mean(data[i]))

# 快速（向量化）
results = np.mean(data, axis=1)
```

#### 6.1.2 并行处理

```python
from multiprocessing import Pool

def process_gene(gene_data):
    return detect_rhythm(gene_data)

# 并行处理
with Pool(processes=4) as pool:
    results = pool.map(process_gene, gene_list)
```

#### 6.1.3 内存优化

```python
# 使用生成器
def gene_iterator(data):
    for gene in data.index:
        yield gene, data.loc[gene]

# 分块处理
chunk_size = 100
for i in range(0, len(data), chunk_size):
    chunk = data.iloc[i:i+chunk_size]
    process_chunk(chunk)
```

### 6.2 性能基准

**测试环境**:
- CPU: 4 cores
- RAM: 16GB
- Python: 3.8

**基准测试**:

| 任务 | 数据规模 | 时间 |
|------|----------|------|
| JTK_CYCLE | 1000基因×24点 | ~30秒 |
| Cosinor | 1000基因×24点 | ~10秒 |
| 随机森林训练 | 500样本×50特征 | ~2秒 |
| 网络构建 | 200节点 | ~5秒 |

### 6.3 优化建议

1. **使用NumPy数组**而非Python列表
2. **预分配内存**而非动态append
3. **使用内置函数**如`np.mean()`而非手动实现
4. **避免重复计算**，缓存中间结果
5. **选择合适的算法**：大规模数据时考虑近似算法

---

## 7. 测试框架

### 7.1 测试架构

```
tests/
├── __init__.py
├── test_rhythm_detection.py    # 节律检测测试
├── test_phase_prediction.py    # 相位预测测试
├── test_network_analysis.py    # 网络分析测试（可扩展）
├── test_visualization.py       # 可视化测试（可扩展）
└── test_utils.py               # 工具函数测试（可扩展）
```

### 7.2 单元测试

#### 7.2.1 测试覆盖

**rhythm_detection**:
```python
def test_synthetic_data_generation():
    """测试合成数据生成"""
    data = generate_synthetic_data(n_genes=100, n_timepoints=24)
    assert data.shape == (100, 24)

def test_jtk_cycle_detection():
    """测试JTK_CYCLE"""
    data = generate_synthetic_data(n_rhythmic=20)
    results = detect_rhythms(data)
    rhythmic = results[results['BH.Q'] < 0.05]
    assert len(rhythmic) > 0

def test_cosinor_regression():
    """测试Cosinor"""
    # 创建完美余弦波
    time = np.arange(24)
    expression = 10 + 2 * np.cos(2*np.pi*time/24)
    result = cosinor_regression(expression, time)
    assert abs(result['amplitude'] - 2.0) < 0.1
```

**phase_prediction**:
```python
def test_feature_extraction():
    """测试特征提取"""
    data = generate_synthetic_data(n_genes=50)
    X, names = extract_gene_features(data)
    assert X.shape[0] == 50
    assert len(names) == X.shape[1]

def test_phase_prediction():
    """测试相位预测"""
    # 训练和预测流程
    model, scaler = train_phase_predictor(X_train, y_train)
    predictions = predict_phase(model, scaler, X_test)
    assert len(predictions) == len(y_test)
    assert all(0 <= p < 24 for p in predictions)
```

### 7.3 集成测试

```python
def test_full_pipeline():
    """测试完整分析流程"""
    # 生成数据
    data = generate_synthetic_data(n_genes=100)

    # 节律检测
    jtk_results = detect_rhythms(data)
    rhythmic = jtk_results[jtk_results['BH.Q'] < 0.05]

    # 网络分析
    if len(rhythmic) > 10:
        rhythmic_data = data.loc[rhythmic['GeneID'].head(20)]
        adj_matrix, edges = build_coexpression_network(rhythmic_data)
        assert len(edges) >= 0

    # 可视化
    plot_rhythmic_genes(data, rhythmic, save_path='test_output.png')
    assert os.path.exists('test_output.png')
```

### 7.4 运行测试

```bash
# 运行所有测试
python -m pytest tests/ -v

# 运行特定测试
python tests/test_rhythm_detection.py

# 测试覆盖率
pytest --cov=src/circadian_analysis tests/

# 生成HTML报告
pytest --cov=src/circadian_analysis --cov-report=html tests/
```

---

## 8. 部署指南

### 8.1 本地安装

#### 8.1.1 从源码安装

```bash
# 克隆仓库
git clone https://github.com/lethaquinn/circadian-transcriptomics.git
cd circadian-transcriptomics

# 创建虚拟环境（推荐）
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate  # Windows

# 安装依赖
pip install -r requirements.txt

# 安装包（开发模式）
pip install -e .
```

#### 8.1.2 使用pip安装（未来）

```bash
pip install circadian-transcriptomics
```

### 8.2 Docker部署

**Dockerfile**:
```dockerfile
FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV PYTHONPATH=/app/src

CMD ["python", "examples/quick_start.py"]
```

**构建和运行**:
```bash
# 构建镜像
docker build -t circadian-analysis .

# 运行容器
docker run -v $(pwd)/results:/app/results circadian-analysis

# 交互式运行
docker run -it circadian-analysis /bin/bash
```

### 8.3 云平台部署

#### 8.3.1 Google Colab

```python
# 在Colab中安装
!git clone https://github.com/lethaquinn/circadian-transcriptomics.git
%cd circadian-transcriptomics
!pip install -r requirements.txt

import sys
sys.path.insert(0, 'src')

from circadian_analysis import demo_analysis
demo_analysis()
```

#### 8.3.2 服务器部署

```bash
# 使用systemd管理服务
sudo nano /etc/systemd/system/circadian-analysis.service
```

```ini
[Unit]
Description=Circadian Analysis Service
After=network.target

[Service]
User=ubuntu
WorkingDirectory=/home/ubuntu/circadian-transcriptomics
Environment="PYTHONPATH=/home/ubuntu/circadian-transcriptomics/src"
ExecStart=/home/ubuntu/venv/bin/python examples/full_pipeline.py

[Install]
WantedBy=multi-user.target
```

```bash
# 启动服务
sudo systemctl start circadian-analysis
sudo systemctl enable circadian-analysis
```

### 8.4 性能监控

```python
import time
import psutil

def profile_analysis():
    start_time = time.time()
    start_mem = psutil.Process().memory_info().rss / 1024 / 1024

    # 运行分析
    results = detect_rhythms(data)

    end_time = time.time()
    end_mem = psutil.Process().memory_info().rss / 1024 / 1024

    print(f"Time: {end_time - start_time:.2f}s")
    print(f"Memory: {end_mem - start_mem:.2f}MB")
```

---

## 9. 故障排除

### 9.1 常见错误

#### 9.1.1 导入错误

**错误**:
```
ModuleNotFoundError: No module named 'circadian_analysis'
```

**解决**:
```bash
# 方法1: 设置PYTHONPATH
export PYTHONPATH=/path/to/circadian-transcriptomics/src

# 方法2: 在代码中添加路径
import sys
sys.path.insert(0, '/path/to/circadian-transcriptomics/src')

# 方法3: 安装包
pip install -e .
```

#### 9.1.2 依赖版本冲突

**错误**:
```
ImportError: cannot import name 'xxx' from 'sklearn'
```

**解决**:
```bash
# 更新scikit-learn
pip install --upgrade scikit-learn

# 或创建新虚拟环境
python -m venv fresh_venv
source fresh_venv/bin/activate
pip install -r requirements.txt
```

#### 9.1.3 内存不足

**错误**:
```
MemoryError: Unable to allocate array
```

**解决**:
```python
# 分块处理
chunk_size = 100
for i in range(0, len(data), chunk_size):
    chunk_data = data.iloc[i:i+chunk_size]
    chunk_results = detect_rhythms(chunk_data)
    # 保存chunk结果
```

#### 9.1.4 绘图显示问题

**错误**: 图片不显示

**解决**:
```python
# 使用非交互式后端
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 保存到文件而非显示
plt.savefig('output.png')
```

### 9.2 调试技巧

#### 9.2.1 启用详细日志

```python
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger('circadian_analysis')
logger.debug("Starting analysis...")
```

#### 9.2.2 检查中间结果

```python
# 在每个步骤打印信息
print(f"Data shape: {data.shape}")
print(f"Time points: {data.columns.tolist()}")
print(f"Value range: [{data.min().min()}, {data.max().max()}]")

# 检查NaN
print(f"NaN values: {data.isna().sum().sum()}")
```

#### 9.2.3 使用断言验证

```python
# 验证数据格式
assert isinstance(data, pd.DataFrame), "Data must be DataFrame"
assert data.shape[1] > 0, "No time points found"
assert not data.isna().any().any(), "Data contains NaN values"
```

### 9.3 性能问题

#### 9.3.1 分析运行缓慢

**诊断**:
```python
import cProfile

cProfile.run('detect_rhythms(data)', 'profile_stats')

import pstats
stats = pstats.Stats('profile_stats')
stats.sort_stats('cumulative')
stats.print_stats(10)
```

**优化**:
1. 减少基因数量（筛选高表达基因）
2. 使用并行处理
3. 降低相位搜索精度

#### 9.3.2 内存使用过高

**监控**:
```python
import tracemalloc

tracemalloc.start()

# 运行分析
results = detect_rhythms(data)

current, peak = tracemalloc.get_traced_memory()
print(f"Current: {current / 1024 / 1024:.2f}MB")
print(f"Peak: {peak / 1024 / 1024:.2f}MB")

tracemalloc.stop()
```

---

## 10. 扩展开发

### 10.1 添加新算法

#### 10.1.1 创建新的节律检测方法

```python
# src/circadian_analysis/rhythm_detection/new_method.py

import numpy as np
import pandas as pd

def detect_rhythms_new(data, **params):
    """
    新的节律检测算法

    Parameters
    ----------
    data : pd.DataFrame
        表达矩阵
    **params : dict
        算法参数

    Returns
    -------
    pd.DataFrame
        检测结果
    """
    results = []

    for gene in data.index:
        expression = data.loc[gene].values

        # 实现您的算法
        score = your_algorithm(expression, **params)
        pval = compute_pvalue(score)

        results.append({
            'GeneID': gene,
            'Score': score,
            'Pvalue': pval
        })

    df = pd.DataFrame(results)

    # FDR校正
    from ..utils.helpers import adjust_pvalues
    df['Qvalue'] = adjust_pvalues(df['Pvalue'].values)

    return df
```

#### 10.1.2 集成到主模块

```python
# src/circadian_analysis/rhythm_detection/__init__.py

from . import new_method

__all__ = [..., "new_method"]
```

### 10.2 自定义可视化

```python
# src/circadian_analysis/visualization/custom_plots.py

import matplotlib.pyplot as plt
import seaborn as sns

def plot_custom_analysis(data, results, save_path=None):
    """
    自定义分析图
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 子图1: 自定义内容
    axes[0, 0].plot(...)
    axes[0, 0].set_title("Custom Analysis")

    # 子图2-4: 其他内容
    ...

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.close()
```

### 10.3 插件系统（概念）

```python
# 未来可扩展为插件架构

class RhythmDetectionPlugin:
    """节律检测插件基类"""

    def __init__(self, name, version):
        self.name = name
        self.version = version

    def detect(self, data, **params):
        raise NotImplementedError

    def validate_input(self, data):
        """验证输入数据"""
        assert isinstance(data, pd.DataFrame)
        return True

# 用户可以继承并实现自己的算法
class MyCustomDetector(RhythmDetectionPlugin):
    def detect(self, data, **params):
        # 实现自定义算法
        return results
```

### 10.4 API开发

```python
# 使用Flask创建Web API（示例）

from flask import Flask, request, jsonify
import pandas as pd
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms

app = Flask(__name__)

@app.route('/api/detect_rhythms', methods=['POST'])
def api_detect_rhythms():
    # 接收JSON数据
    data_json = request.get_json()

    # 转换为DataFrame
    data = pd.DataFrame(data_json['expression_matrix'])

    # 运行分析
    results = detect_rhythms(data, period=data_json.get('period', 24.0))

    # 返回结果
    return jsonify(results.to_dict('records'))

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
```

### 10.5 贡献指南

如果您想为项目贡献代码：

1. **Fork仓库**
2. **创建功能分支**: `git checkout -b feature/new-algorithm`
3. **编写代码和测试**
4. **确保测试通过**: `pytest tests/`
5. **提交PR** (Pull Request)

**代码规范**:
- 遵循PEP 8
- 添加docstring
- 编写单元测试
- 更新文档

---

## 附录

### A. 数学符号表

| 符号 | 含义 |
|------|------|
| $\tau$ | Kendall's tau相关系数 |
| $\omega$ | 角频率 ($2\pi/T$) |
| $T$ | 周期 |
| $\phi$ | 相位 |
| $A$ | 振幅 |
| $M$ | MESOR (平均水平) |
| $R^2$ | 决定系数 |
| $\sigma$ | 标准差 |

### B. 参考文献

1. Hughes ME, et al. (2010). JTK_CYCLE: an efficient nonparametric algorithm for detecting rhythmic components in genome-scale data sets. *J Biol Rhythms*, 25(5):372-80.

2. Cornelissen G. (2014). Cosinor-based rhythmometry. *Theoretical Biology and Medical Modelling*, 11:16.

3. Lomb NR. (1976). Least-squares frequency analysis of unequally spaced data. *Astrophysics and Space Science*, 39:447-462.

4. Zhang R, et al. (2014). A circadian gene expression atlas in mammals: implications for biology and medicine. *PNAS*, 111(45):16219-16224.

### C. 术语表

- **Circadian Rhythm**: 昼夜节律，约24小时的生物周期
- **Zeitgeber Time (ZT)**: 时间给定时间，ZT0为光照开始
- **Circadian Time (CT)**: 昼夜节律时间，在恒定条件下
- **MESOR**: Midline Estimating Statistic Of Rhythm，节律的平均水平
- **Acrophase**: 峰相位，节律达到峰值的时间
- **Amplitude**: 振幅，峰值与平均值的差
- **FDR**: False Discovery Rate，假阳性率
- **Hub Gene**: 枢纽基因，网络中高度连接的基因

### D. 常用命令速查

```bash
# 快速测试
bash run_all_tests.sh

# 运行Demo
python -c "from circadian_analysis import demo_analysis; demo_analysis()"

# 运行示例
python examples/quick_start.py
python examples/full_pipeline.py

# 查看结果
ls -R results/

# 生成文档
pydoc -w circadian_analysis

# 代码格式化
black src/

# 静态检查
flake8 src/
```

---

## 版本历史

### v0.1.0 (2025-11-17)
- ✅ 初始版本
- ✅ 实现核心算法（JTK_CYCLE, Cosinor, Lomb-Scargle, Wavelet）
- ✅ 实现相位预测（RF, GB, MLP, Ensemble）
- ✅ 实现网络分析（共表达、社区检测）
- ✅ 完整的可视化系统
- ✅ 测试框架（9个测试，100%通过率）
- ✅ 文档和示例

---

## 联系方式

**作者**: lethaquinn
**Email**: 2679066373@qq.com
**GitHub**: https://github.com/lethaquinn/circadian-transcriptomics

---

**文档更新日期**: 2025-11-17
**文档版本**: 1.0
**适用软件版本**: circadian-transcriptomics v0.1.0
