# 📖 用户指南

**Circadian Transcriptomics Analysis Pipeline - 使用手册**

欢迎使用昼夜节律转录组学分析管道！本指南将帮助您快速上手并充分利用这个工具包。

---

## 🎯 快速导航

- [安装](#安装)
- [快速开始](#快速开始)
- [数据准备](#数据准备)
- [工作流程](#工作流程)
- [常见应用场景](#常见应用场景)
- [最佳实践](#最佳实践)
- [FAQ](#faq)

---

## 📦 安装

### 系统要求

- Python ≥ 3.8
- 8GB RAM（推荐16GB用于大规模数据）
- 2GB 磁盘空间

### 安装步骤

```bash
# 1. 克隆仓库
git clone https://github.com/lethaquinn/circadian-transcriptomics.git
cd circadian-transcriptomics

# 2. 安装依赖
pip install -r requirements.txt

# 3. 设置Python路径
export PYTHONPATH=$PWD/src

# 4. 测试安装
python -c "from circadian_analysis import demo_analysis; demo_analysis()"
```

看到"✅ Demo completed successfully!"说明安装成功！

---

## 🚀 快速开始

### 30秒入门

```python
# 导入包
import sys
sys.path.insert(0, 'src')

from circadian_analysis import demo_analysis

# 运行演示
demo_analysis()
```

这会：
1. 生成合成昼夜节律数据
2. 检测节律基因
3. 创建可视化图表
4. 保存结果到 `results/demo/`

### 5分钟完整示例

```python
from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.visualization.phase_plots import (
    plot_rhythmic_genes, plot_phase_distribution
)

# 1. 准备数据（这里使用合成数据，实际使用时替换为您的数据）
data = generate_synthetic_data(n_genes=200, n_timepoints=24, n_rhythmic=50)

# 2. 检测节律
results = detect_rhythms(data, period=24.0)

# 3. 筛选显著基因
rhythmic_genes = results[results['BH.Q'] < 0.05]
print(f"检测到 {len(rhythmic_genes)} 个节律基因")

# 4. 可视化
plot_rhythmic_genes(data, rhythmic_genes, save_path='my_results.png')
plot_phase_distribution(rhythmic_genes, save_path='phase_dist.png')

# 5. 保存结果
results.to_csv('rhythm_results.csv', index=False)
```

---

## 📊 数据准备

### 数据格式要求

您的数据应该是一个**基因 × 时间点**的矩阵：

```
           ZT00   ZT02   ZT04   ZT06  ...  ZT22
Gene_001   5.23   5.67   6.12   7.89  ...  4.98
Gene_002  10.50  10.20   9.87   9.45  ...  10.80
Gene_003   3.45   3.12   2.98   3.34  ...  3.67
...
```

### 从文件加载数据

#### CSV文件

```python
import pandas as pd

# 读取CSV
data = pd.read_csv('expression_data.csv', index_col=0)

# 检查格式
print(f"数据维度: {data.shape}")  # (基因数, 时间点数)
print(f"列名: {data.columns.tolist()}")
```

#### Excel文件

```python
data = pd.read_excel('expression_data.xlsx', index_col=0, sheet_name='Sheet1')
```

### 从GEO下载数据

```python
from data.download_scripts.geo_downloader import download_geo_dataset

# 下载GEO数据集
file_path = download_geo_dataset('GSE11923', output_dir='data/raw')

# 解析并加载
# (需要进一步处理，具体取决于数据格式)
```

### 数据预处理

```python
from data.download_scripts.data_preprocessor import preprocess_pipeline

# 完整预处理流程
processed_data = preprocess_pipeline(
    raw_data,
    filter_threshold=1.0,    # 过滤低表达基因
    normalize=True,          # 分位数标准化
    log_transform=True       # Log2转换
)
```

---

## 🔬 工作流程

### 流程1: 基础节律检测

**目标**: 找出哪些基因具有昼夜节律

```python
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms

# Step 1: 运行JTK_CYCLE
results = detect_rhythms(data, period=24.0, alpha=0.05)

# Step 2: 筛选显著基因
rhythmic = results[results['BH.Q'] < 0.05]

# Step 3: 查看结果
print(f"节律基因数: {len(rhythmic)}")
print(f"节律基因比例: {100*len(rhythmic)/len(results):.1f}%")

# Step 4: 导出节律基因列表
rhythmic['GeneID'].to_csv('rhythmic_genes.txt', index=False, header=False)
```

**结果解读**:
- `BH.Q < 0.05`: 该基因具有显著昼夜节律
- `LAG`: 该基因表达峰值的时间（小时）
- `AMP`: 振幅，表达变化的幅度
- `TAU`: Kendall相关系数，值越大节律越强

### 流程2: 比较多种方法

**目标**: 使用多种算法交叉验证

```python
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.rhythm_detection.cosinor import detect_rhythms_cosinor
from circadian_analysis.visualization.dashboard import create_comparison_dashboard

# 运行多种方法
jtk_results = detect_rhythms(data)
cosinor_results = detect_rhythms_cosinor(data)

# 比较结果
results_dict = {
    'JTK_CYCLE': jtk_results,
    'Cosinor': cosinor_results
}

# 创建对比图
create_comparison_dashboard(results_dict, save_path='comparison.png')

# 找出两种方法都检测到的基因（高可信度）
jtk_genes = set(jtk_results[jtk_results['BH.Q'] < 0.05]['GeneID'])
cosinor_genes = set(cosinor_results[cosinor_results['Qvalue'] < 0.05]['GeneID'])

high_confidence = jtk_genes & cosinor_genes
print(f"高可信度节律基因: {len(high_confidence)}")
```

### 流程3: 相位分析

**目标**: 分析基因表达的峰值时间

```python
from circadian_analysis.visualization.phase_plots import (
    plot_phase_distribution, plot_phase_wheel
)

# 相位分布
plot_phase_distribution(
    rhythmic_results,
    phase_col='LAG',
    save_path='phase_distribution.png'
)

# 相位轮盘图
plot_phase_wheel(
    rhythmic_results,
    top_n=50,
    save_path='phase_wheel.png'
)

# 按相位分组
import numpy as np

def assign_phase_group(phase):
    """将相位分配到时间段"""
    if 0 <= phase < 6:
        return '凌晨 (0-6h)'
    elif 6 <= phase < 12:
        return '上午 (6-12h)'
    elif 12 <= phase < 18:
        return '下午 (12-18h)'
    else:
        return '晚上 (18-24h)'

rhythmic_results['PhaseGroup'] = rhythmic_results['LAG'].apply(assign_phase_group)

# 统计各时间段基因数
phase_counts = rhythmic_results['PhaseGroup'].value_counts()
print(phase_counts)
```

### 流程4: 网络分析

**目标**: 找出相互作用的节律基因

```python
from circadian_analysis.network_analysis.coexpression import (
    build_coexpression_network, identify_hub_genes
)
from circadian_analysis.network_analysis.community import (
    greedy_modularity_communities, annotate_communities
)

# 只分析节律基因
rhythmic_data = data.loc[rhythmic_results['GeneID'].values]

# 构建共表达网络
adj_matrix, edge_list = build_coexpression_network(
    rhythmic_data,
    threshold=0.8,  # 只保留高相关的连接
    method='pearson'
)

print(f"网络: {len(adj_matrix)} 节点, {len(edge_list)} 边")

# 找出Hub基因（高度连接的关键基因）
hubs = identify_hub_genes(adj_matrix, top_n=10)
print("\nTop 10 Hub基因:")
print(hubs[['Gene', 'Degree']])

# 社区检测（功能模块）
communities = greedy_modularity_communities(adj_matrix)
comm_df = annotate_communities(communities)

print(f"\n检测到 {comm_df['Community'].nunique()} 个功能模块")

# 导出结果
edge_list.to_csv('network_edges.csv', index=False)
hubs.to_csv('hub_genes.csv', index=False)
comm_df.to_csv('communities.csv', index=False)
```

### 流程5: 相位预测

**目标**: 基于表达模式预测未知基因的相位

```python
from circadian_analysis.phase_prediction.ml_models import (
    extract_gene_features, train_phase_predictor, predict_phase
)
from sklearn.model_selection import train_test_split

# 准备训练数据（使用已知节律的基因）
known_rhythmic = data.loc[rhythmic_results['GeneID'].head(100).values]
X, feature_names = extract_gene_features(known_rhythmic, rhythmic_results.head(100))
y = rhythmic_results.head(100)['LAG'].values

# 划分训练/测试集
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# 训练模型
model, scaler = train_phase_predictor(
    X_train, y_train,
    model_type='random_forest',
    n_estimators=200
)

# 在测试集上评估
y_pred = predict_phase(model, scaler, X_test)

from circadian_analysis.utils.metrics import compute_phase_error
error = compute_phase_error(y_test, y_pred)
print(f"预测误差: {error:.2f} 小时")

# 预测未知基因的相位
unknown_genes = data.loc[other_gene_ids]
X_unknown, _ = extract_gene_features(unknown_genes)
predicted_phases = predict_phase(model, scaler, X_unknown)

# 保存预测结果
pd.DataFrame({
    'GeneID': other_gene_ids,
    'PredictedPhase': predicted_phases
}).to_csv('predicted_phases.csv', index=False)
```

---

## 💡 常见应用场景

### 场景1: 发现新的昼夜节律基因

**问题**: 我有一个转录组数据集，想知道哪些基因具有昼夜节律

**解决方案**:
```python
# 1. 加载数据
data = pd.read_csv('my_data.csv', index_col=0)

# 2. 运行检测
results = detect_rhythms(data, period=24.0)

# 3. 筛选并排序
rhythmic = results[results['BH.Q'] < 0.05].sort_values('BH.Q')

# 4. 保存Top100节律基因
top100 = rhythmic.head(100)
top100.to_csv('top100_rhythmic_genes.csv')

# 5. 可视化
plot_rhythmic_genes(data, top100, top_n=12, save_path='top12.png')
```

### 场景2: 比较不同条件下的节律

**问题**: 我有对照组和处理组的数据，想比较节律的变化

**解决方案**:
```python
# 分别分析两组
control_results = detect_rhythms(control_data)
treatment_results = detect_rhythms(treatment_data)

# 找出节律改变的基因
control_rhythmic = set(control_results[control_results['BH.Q'] < 0.05]['GeneID'])
treatment_rhythmic = set(treatment_results[treatment_results['BH.Q'] < 0.05]['GeneID'])

# 新获得节律的基因
gained_rhythm = treatment_rhythmic - control_rhythmic

# 失去节律的基因
lost_rhythm = control_rhythmic - treatment_rhythmic

print(f"获得节律: {len(gained_rhythm)} 个基因")
print(f"失去节律: {len(lost_rhythm)} 个基因")

# 相位偏移分析
merged = control_results.merge(
    treatment_results,
    on='GeneID',
    suffixes=('_Control', '_Treatment')
)

# 只看两组都显著的基因
both_sig = merged[
    (merged['BH.Q_Control'] < 0.05) &
    (merged['BH.Q_Treatment'] < 0.05)
]

# 计算相位偏移
from circadian_analysis.network_analysis.phase_coupling import compute_phase_difference
phase_shift = compute_phase_difference(
    both_sig['LAG_Treatment'].values,
    both_sig['LAG_Control'].values
)

both_sig['PhaseShift'] = phase_shift

# 找出相位显著偏移的基因（>3小时）
phase_shifted = both_sig[abs(both_sig['PhaseShift']) > 3]
print(f"相位显著偏移: {len(phase_shifted)} 个基因")
```

### 场景3: 功能富集分析

**问题**: 节律基因富集在哪些生物学功能？

**解决方案**:
```python
# 导出节律基因列表
rhythmic_genes = rhythmic_results['GeneID'].tolist()

# 按相位分组导出
for phase_group in ['凌晨', '上午', '下午', '晚上']:
    genes = rhythmic_results[
        rhythmic_results['PhaseGroup'] == f'{phase_group} (X-Yh)'
    ]['GeneID'].tolist()

    with open(f'genes_{phase_group}.txt', 'w') as f:
        f.write('\n'.join(genes))

# 然后使用这些基因列表进行GO/KEGG富集分析
# 可以使用在线工具如DAVID, Metascape, g:Profiler等
```

### 场景4: 时间序列聚类

**问题**: 将节律基因按表达模式分组

**解决方案**:
```python
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

# 提取节律基因的表达数据
rhythmic_data = data.loc[rhythmic_results['GeneID'].values]

# 标准化
from circadian_analysis.utils.helpers import normalize_expression
normalized = normalize_expression(rhythmic_data, method='zscore')

# K-means聚类
n_clusters = 6
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
clusters = kmeans.fit_predict(normalized.values)

# 添加聚类标签
rhythmic_results['Cluster'] = clusters

# 可视化每个cluster的平均表达模式
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for i in range(n_clusters):
    cluster_data = normalized[clusters == i]
    mean_pattern = cluster_data.mean(axis=0)

    axes[i].plot(range(len(mean_pattern)), mean_pattern, linewidth=2)
    axes[i].set_title(f'Cluster {i+1} (n={sum(clusters==i)})')
    axes[i].set_xlabel('时间点')
    axes[i].set_ylabel('Z-score')
    axes[i].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('expression_clusters.png', dpi=300)

# 保存每个cluster的基因列表
for i in range(n_clusters):
    cluster_genes = rhythmic_results[rhythmic_results['Cluster'] == i]['GeneID']
    cluster_genes.to_csv(f'cluster_{i+1}_genes.txt', index=False, header=False)
```

---

## ✨ 最佳实践

### 1. 数据质量控制

**在分析前检查**:
```python
# 检查缺失值
print(f"缺失值: {data.isna().sum().sum()}")

# 检查数据范围
print(f"表达范围: [{data.min().min():.2f}, {data.max().max():.2f}]")

# 检查分布
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
data.values.flatten().hist(bins=50)
plt.xlabel('表达值')
plt.ylabel('频数')
plt.title('表达值分布')
plt.show()
```

### 2. 选择合适的算法

| 场景 | 推荐算法 | 原因 |
|------|----------|------|
| 标准24h数据 | JTK_CYCLE | 非参数，鲁棒性强 |
| 需要量化振幅和相位 | Cosinor | 提供精确的参数估计 |
| 非等间隔采样 | Lomb-Scargle | 专门处理不规则采样 |
| 周期未知 | Lomb-Scargle / Wavelet | 可搜索多个周期 |
| 周期变化的数据 | Wavelet | 时频联合分析 |

### 3. 参数选择

```python
# JTK_CYCLE
results = detect_rhythms(
    data,
    period=24.0,        # 如果是小鼠数据可能需要调整
    lag_range=(0, 23),  # 测试所有可能相位
    alpha=0.05          # FDR阈值，可以调整为0.01更严格
)

# Cosinor
results = detect_rhythms_cosinor(
    data,
    period=24.0,  # 也可以尝试[23.0, 24.0, 25.0]范围
)

# 共表达网络
adj_matrix, edges = build_coexpression_network(
    data,
    threshold=0.7,      # 阈值越高，网络越稀疏但连接更可靠
    method='pearson'    # 或'spearman'用于非线性关系
)
```

### 4. 可视化技巧

```python
# 高质量发表图
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'

# 保存为矢量图（可编辑）
plt.savefig('figure.pdf', format='pdf', bbox_inches='tight')

# 或高分辨率位图
plt.savefig('figure.png', dpi=600, bbox_inches='tight')
```

### 5. 结果验证

**交叉验证**:
```python
# 使用多种方法
jtk_genes = set(jtk_results[jtk_results['BH.Q'] < 0.05]['GeneID'])
cosinor_genes = set(cosinor_results[cosinor_results['Qvalue'] < 0.05]['GeneID'])

# 计算一致性
overlap = len(jtk_genes & cosinor_genes)
union = len(jtk_genes | cosinor_genes)
jaccard = overlap / union

print(f"Jaccard相似度: {jaccard:.2f}")
```

**文献验证**:
```python
# 检查已知昼夜节律基因是否被检测到
known_clock_genes = ['Per1', 'Per2', 'Cry1', 'Cry2', 'Bmal1', 'Clock']

detected_clock = [g for g in known_clock_genes if g in rhythmic_results['GeneID'].values]
print(f"检测到的时钟基因: {detected_clock}")
```

---

## ❓ FAQ

### Q1: 我的数据只有12个时间点，可以用吗？

**A**: 可以，但建议至少有**12个时间点**覆盖一个完整周期。数据点越多，检测越准确。如果只有少量时间点：
- JTK_CYCLE仍然可用，但功效降低
- Cosinor更稳定，推荐使用
- 建议降低显著性阈值（如alpha=0.1）

### Q2: 我的数据是counts，需要标准化吗？

**A**: 是的，强烈建议标准化：

```python
import numpy as np

# Log2转换
data_log = np.log2(data + 1)

# 或使用预处理pipeline
from data.download_scripts.data_preprocessor import preprocess_pipeline
data_processed = preprocess_pipeline(data, log_transform=True, normalize=True)
```

### Q3: 检测到的节律基因太少怎么办？

**可能原因和解决方案**:

1. **数据质量问题**
   - 检查是否标准化
   - 检查是否有异常值
   - 尝试不同的标准化方法

2. **参数过于严格**
   - 放宽FDR阈值：`alpha=0.1`
   - 尝试不同周期：`period=[23.0, 24.0, 25.0]`

3. **实际确实节律基因少**
   - 某些组织/条件下节律基因本就较少
   - 这是正常的生物学现象

### Q4: 如何处理重复样本？

**A**: 有两种方式：

```python
# 方法1: 平均重复样本
averaged_data = replicate_data.groupby(level=0, axis=1).mean()

# 方法2: 使用预处理工具
from data.download_scripts.data_preprocessor import handle_replicates

processed = handle_replicates(
    data,
    replicate_info,  # DataFrame标注哪些样本是重复
    method='mean'    # 或'median'
)
```

### Q5: 结果不同算法差异很大？

**A**: 这很正常，不同算法有不同假设：

- **JTK_CYCLE**: 非参数，检测任何单调趋势
- **Cosinor**: 假设余弦波形
- **Lomb-Scargle**: 类似Cosinor，但更灵活

**建议**:
1. 关注多种方法都检测到的高可信度基因
2. 根据生物学背景选择合适方法
3. 报告时说明使用的方法和参数

### Q6: 可以分析RNA-seq数据吗？

**A**: 可以！RNA-seq是转录组学的标准方法。

**注意事项**:
```python
# RNA-seq counts需要标准化
# 推荐TPM或FPKM，或使用DESeq2/edgeR标准化

# 示例：简单标准化
# 1. 过滤低表达基因
data_filtered = data[data.sum(axis=1) > 10]

# 2. CPM标准化
total_counts = data_filtered.sum(axis=0)
data_cpm = (data_filtered / total_counts) * 1e6

# 3. Log转换
data_log = np.log2(data_cpm + 1)

# 4. Z-score标准化（per gene）
data_zscore = (data_log.T - data_log.mean(axis=1)) / data_log.std(axis=1)
data_zscore = data_zscore.T
```

### Q7: 如何解释相位(Phase/LAG)?

**A**: 相位表示基因表达达到峰值的时间：

- `LAG = 0` → 0点（午夜）达到峰值
- `LAG = 6` → 6点（清晨）达到峰值
- `LAG = 12` → 12点（中午）达到峰值
- `LAG = 18` → 18点（傍晚）达到峰值

**ZT (Zeitgeber Time)**:
- ZT0 = 光照开始（通常是"dawn"）
- ZT12 = 黑暗开始（通常是"dusk"）

### Q8: 网络分析的阈值如何选择？

**A**: 相关系数阈值的选择：

```python
# 方法1: 固定阈值
threshold = 0.8  # 高阈值 = 稀疏但高质量的网络

# 方法2: Top-k边
import numpy as np
corr_matrix = data.corr()
threshold = np.percentile(corr_matrix.values.flatten(), 95)  # Top 5%

# 方法3: 根据网络密度
# 尝试不同阈值，选择合适的网络密度（如边数/最大可能边数 ≈ 0.01-0.05）
```

**经验法则**:
- 0.5-0.6: 很宽松，网络密集
- 0.7-0.8: 适中
- 0.9+: 很严格，只保留最强关联

---

## 📚 进一步学习

- **技术文档**: `docs/TECHNICAL_DOCUMENTATION.md` - 深入的技术细节
- **API参考**: `docs/API_REFERENCE.md` - 所有函数的详细说明
- **示例代码**: `examples/` - 更多应用示例
- **测试代码**: `tests/` - 看看我们如何测试每个功能

---

## 🆘 获取帮助

遇到问题？

1. 查看 [FAQ](#faq)
2. 检查 [故障排除](TECHNICAL_DOCUMENTATION.md#故障排除)
3. 查看示例代码 `examples/`
4. 提issue到GitHub仓库

---

## 🎓 引用

如果您在研究中使用了本工具，请引用：

```
Lethaquinn (2025). Circadian Transcriptomics Analysis Pipeline:
Advanced computational framework for circadian rhythm and biomarker discovery.
GitHub: https://github.com/lethaquinn/circadian-transcriptomics
```

---

**祝您分析顺利！**  🚀

如果本工具对您有帮助，请给我们一个⭐Star！
