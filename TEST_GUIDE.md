# 🧪 Circadian Transcriptomics - 測試指南

本指南將幫助您驗證項目的所有功能是否正常運行。

## 📋 目錄
1. [環境檢查](#環境檢查)
2. [快速測試](#快速測試)
3. [單元測試](#單元測試)
4. [功能測試](#功能測試)
5. [完整測試](#完整測試)
6. [進階測試](#進階測試)

---

## 1️⃣ 環境檢查

### 檢查Python版本
```bash
python --version  # 需要 Python 3.8+
```

### 安裝依賴
```bash
pip install -r requirements.txt
```

### 驗證安裝
```bash
python -c "import numpy, pandas, scipy, sklearn, matplotlib; print('✅ 所有依賴已安裝')"
```

---

## 2️⃣ 快速測試（5分鐘）

### 測試 1: 最簡單的導入測試
```bash
cd /home/user/circadian-transcriptomics
export PYTHONPATH=src
python -c "from circadian_analysis import demo_analysis; print('✅ 模組導入成功')"
```

### 測試 2: 運行內建Demo
```bash
python -c "from circadian_analysis import demo_analysis; demo_analysis()"
```

**預期結果**:
- 生成合成數據
- 檢測節律基因
- 保存圖表到 `results/demo/`
- 顯示 "✅ Demo completed successfully!"

### 測試 3: 驗證輸出文件
```bash
ls -lh results/demo/
# 應該看到: rhythmic_genes.png
```

---

## 3️⃣ 單元測試（10分鐘）

### 運行所有測試
```bash
# 測試節律檢測
python tests/test_rhythm_detection.py

# 測試相位預測
python tests/test_phase_prediction.py
```

**預期結果**: 所有測試顯示 ✓

### 使用pytest運行（推薦）
```bash
pip install pytest
pytest tests/ -v
```

---

## 4️⃣ 功能測試（15分鐘）

### 測試 1: 快速開始示例
```bash
python examples/quick_start.py
```

**檢查點**:
- [ ] 生成200個基因的數據
- [ ] 檢測到節律基因（應該約50-60個）
- [ ] 生成3個輸出文件在 `results/quick_start/`:
  - `top_rhythmic_genes.png` - 展示前6個節律基因
  - `phase_distribution.png` - 相位分佈圖
  - `jtk_results.csv` - 完整結果表

**驗證輸出**:
```bash
ls results/quick_start/
wc -l results/quick_start/jtk_results.csv  # 應該有201行（包含表頭）
```

### 測試 2: 完整流程示例
```bash
python examples/full_pipeline.py
```

**檢查點**:
- [ ] 生成500個基因數據
- [ ] 運行JTK_CYCLE和Cosinor兩種方法
- [ ] 構建共表達網絡
- [ ] 生成綜合儀表板
- [ ] 輸出5個文件在 `results/full_pipeline/`

**驗證輸出**:
```bash
ls results/full_pipeline/
# 應該看到:
# - analysis_dashboard.png
# - jtk_results.csv
# - cosinor_results.csv
# - network_edges.csv
# - network_metrics.csv
```

---

## 5️⃣ 完整測試（30分鐘）

### 測試各個模組功能

#### A. 節律檢測
```bash
cat > test_rhythm.py << 'EOF'
import sys
sys.path.insert(0, 'src')

from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.rhythm_detection.cosinor import detect_rhythms_cosinor

# 生成測試數據
data = generate_synthetic_data(n_genes=100, n_timepoints=24, n_rhythmic=30)

# 測試JTK_CYCLE
print("Testing JTK_CYCLE...")
jtk_results = detect_rhythms(data, period=24.0)
jtk_sig = jtk_results[jtk_results['BH.Q'] < 0.05]
print(f"✅ JTK_CYCLE: 檢測到 {len(jtk_sig)} 個節律基因")

# 測試Cosinor
print("\nTesting Cosinor...")
cosinor_results = detect_rhythms_cosinor(data, period=24.0)
cosinor_sig = cosinor_results[cosinor_results['Qvalue'] < 0.05]
print(f"✅ Cosinor: 檢測到 {len(cosinor_sig)} 個節律基因")

print("\n🎉 節律檢測測試完成!")
EOF

python test_rhythm.py
```

#### B. 相位預測
```bash
cat > test_phase.py << 'EOF'
import sys
sys.path.insert(0, 'src')

from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.phase_prediction.ml_models import (
    extract_gene_features, train_phase_predictor, predict_phase
)
import numpy as np

# 生成數據
data = generate_synthetic_data(n_genes=100, n_timepoints=24, n_rhythmic=50)

# 獲取節律基因
results = detect_rhythms(data, period=24.0)
rhythmic = results[results['BH.Q'] < 0.05].head(50)

# 提取特徵
rhythmic_data = data.loc[rhythmic['GeneID'].values]
X, _ = extract_gene_features(rhythmic_data)
y = rhythmic['LAG'].values

print(f"特徵矩陣: {X.shape}")
print(f"相位標籤: {y.shape}")

# 訓練模型
n_train = int(0.8 * len(X))
X_train, X_test = X[:n_train], X[n_train:]
y_train, y_test = y[:n_train], y[n_train:]

print("\n訓練隨機森林模型...")
model, scaler = train_phase_predictor(X_train, y_train, model_type='random_forest')

# 預測
predictions = predict_phase(model, scaler, X_test)
mae = np.mean(np.abs(predictions - y_test))

print(f"✅ 預測平均誤差: {mae:.2f} 小時")
print("🎉 相位預測測試完成!")
EOF

python test_phase.py
```

#### C. 網絡分析
```bash
cat > test_network.py << 'EOF'
import sys
sys.path.insert(0, 'src')

from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.network_analysis.coexpression import (
    build_coexpression_network, compute_network_metrics
)

# 生成數據
data = generate_synthetic_data(n_genes=50, n_timepoints=24, n_rhythmic=30)

# 獲取節律基因
results = detect_rhythms(data, period=24.0)
rhythmic = results[results['BH.Q'] < 0.05]
rhythmic_data = data.loc[rhythmic['GeneID'].values]

print(f"節律基因數: {len(rhythmic_data)}")

# 構建網絡
print("\n構建共表達網絡...")
adj_matrix, edge_list = build_coexpression_network(
    rhythmic_data, threshold=0.7, method='pearson'
)

print(f"✅ 網絡節點數: {len(adj_matrix)}")
print(f"✅ 網絡邊數: {len(edge_list)}")

if len(edge_list) > 0:
    # 計算網絡指標
    metrics = compute_network_metrics(adj_matrix)
    print(f"✅ 平均度: {metrics['Degree'].mean():.2f}")

print("🎉 網絡分析測試完成!")
EOF

python test_network.py
```

#### D. 可視化
```bash
cat > test_viz.py << 'EOF'
import sys
sys.path.insert(0, 'src')
import os

from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.visualization.phase_plots import (
    plot_rhythmic_genes, plot_phase_distribution
)

# 生成數據
data = generate_synthetic_data(n_genes=100, n_timepoints=24, n_rhythmic=30)

# 檢測節律
results = detect_rhythms(data, period=24.0)
rhythmic = results[results['BH.Q'] < 0.05]

print(f"節律基因數: {len(rhythmic)}")

# 創建輸出目錄
os.makedirs('results/test_viz', exist_ok=True)

# 繪圖
print("\n生成可視化圖表...")
plot_rhythmic_genes(data, rhythmic, top_n=6,
                   save_path='results/test_viz/rhythmic_genes.png')

plot_phase_distribution(rhythmic, phase_col='LAG',
                       save_path='results/test_viz/phase_dist.png')

print("✅ 圖表已保存到 results/test_viz/")
print("🎉 可視化測試完成!")
EOF

python test_viz.py
```

---

## 6️⃣ 進階測試（自定義數據）

### 測試自己的數據

創建一個測試腳本來分析您自己的數據：

```bash
cat > analyze_my_data.py << 'EOF'
import sys
sys.path.insert(0, 'src')
import pandas as pd

from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.visualization.phase_plots import plot_rhythmic_genes

# 讀取您的數據
# 數據格式: 行=基因, 列=時間點
# data = pd.read_csv('your_data.csv', index_col=0)

# 或使用示例數據
from circadian_analysis.utils.helpers import generate_synthetic_data
data = generate_synthetic_data(n_genes=200, n_timepoints=24, n_rhythmic=50)

print(f"數據維度: {data.shape}")

# 檢測節律
print("\n運行節律檢測...")
results = detect_rhythms(data, period=24.0, alpha=0.05)

# 篩選顯著基因
rhythmic = results[results['BH.Q'] < 0.05]
print(f"檢測到 {len(rhythmic)} 個節律基因 (q < 0.05)")

# 保存結果
results.to_csv('my_results.csv', index=False)
print("\n結果已保存到 my_results.csv")

# 可視化
plot_rhythmic_genes(data, rhythmic, top_n=6, save_path='my_plot.png')
print("圖表已保存到 my_plot.png")
EOF

python analyze_my_data.py
```

---

## 📊 測試檢查清單

完成以下所有檢查項即表示測試成功：

### 基礎測試
- [ ] Python環境正確（≥3.8）
- [ ] 所有依賴安裝成功
- [ ] 模組可以正常導入
- [ ] Demo運行成功

### 單元測試
- [ ] test_rhythm_detection.py 全部通過（6個測試）
- [ ] test_phase_prediction.py 全部通過（3個測試）

### 功能測試
- [ ] quick_start.py 運行成功
- [ ] full_pipeline.py 運行成功
- [ ] 生成的CSV文件可讀且格式正確
- [ ] 生成的PNG圖片可查看

### 模組測試
- [ ] JTK_CYCLE檢測功能正常
- [ ] Cosinor檢測功能正常
- [ ] 相位預測模型可訓練
- [ ] 網絡構建功能正常
- [ ] 可視化圖表生成正常

---

## 🐛 常見問題

### 問題1: 導入錯誤
```
ModuleNotFoundError: No module named 'circadian_analysis'
```
**解決方案**:
```bash
export PYTHONPATH=src
# 或在Python中:
import sys
sys.path.insert(0, 'src')
```

### 問題2: 缺少依賴
```
ModuleNotFoundError: No module named 'xxx'
```
**解決方案**:
```bash
pip install -r requirements.txt
```

### 問題3: 圖表不顯示
**解決方案**: 圖表已保存為文件，查看 `results/` 目錄下的PNG文件

### 問題4: 測試運行緩慢
這是正常的，節律檢測需要一些計算時間，特別是基因數量較多時。

---

## 📞 獲取幫助

如果測試過程中遇到問題：

1. 檢查錯誤訊息
2. 確認所有依賴已安裝
3. 查看生成的日誌文件
4. 檢查 `results/` 目錄是否有輸出

---

## ✅ 測試成功標準

完成所有測試後，您應該看到：

- ✅ 9個單元測試全部通過
- ✅ 3個示例程序成功運行
- ✅ `results/` 目錄下有多個子目錄和輸出文件
- ✅ CSV文件包含完整的分析結果
- ✅ PNG圖片展示清晰的可視化結果

**恭喜！您的項目已經完全可以使用了！** 🎉
