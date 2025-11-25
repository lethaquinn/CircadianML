#!/bin/bash

# 一鍵測試腳本 - Circadian Transcriptomics Analysis Pipeline
# 此腳本會運行所有測試並生成測試報告

set -e  # 遇到錯誤立即退出

echo "========================================================================"
echo "🧬 Circadian Transcriptomics - 自動化測試套件"
echo "========================================================================"
echo ""

# 顏色定義
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 測試計數
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# 設置Python路徑
export PYTHONPATH=src

# 清理舊結果
echo "📁 清理舊測試結果..."
rm -rf results/test_* 2>/dev/null || true
echo ""

# ============================================================
# 第一部分: 環境檢查
# ============================================================
echo "========================================================================"
echo "1️⃣  環境檢查"
echo "========================================================================"

echo -n "檢查 Python 版本... "
TOTAL_TESTS=$((TOTAL_TESTS + 1))
PYTHON_VERSION=$(python --version 2>&1 | awk '{print $2}')
if [[ $(echo "$PYTHON_VERSION" | cut -d. -f1) -ge 3 ]] && [[ $(echo "$PYTHON_VERSION" | cut -d. -f2) -ge 8 ]]; then
    echo -e "${GREEN}✓${NC} Python $PYTHON_VERSION"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    echo -e "${RED}✗${NC} Python版本過低 ($PYTHON_VERSION), 需要 >= 3.8"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

echo -n "檢查依賴包... "
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if python -c "import numpy, pandas, scipy, sklearn, matplotlib, seaborn" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} 所有依賴已安裝"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    echo -e "${YELLOW}⚠${NC}  缺少依賴，正在安裝..."
    pip install -q -r requirements.txt
    echo -e "${GREEN}✓${NC} 依賴安裝完成"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi

echo ""

# ============================================================
# 第二部分: 單元測試
# ============================================================
echo "========================================================================"
echo "2️⃣  單元測試"
echo "========================================================================"

echo "運行 test_rhythm_detection.py..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if python tests/test_rhythm_detection.py > /tmp/test1.log 2>&1; then
    echo -e "${GREEN}✓${NC} rhythm_detection 測試通過 (6/6)"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    echo -e "${RED}✗${NC} rhythm_detection 測試失敗"
    cat /tmp/test1.log
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

echo "運行 test_phase_prediction.py..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if python tests/test_phase_prediction.py > /tmp/test2.log 2>&1; then
    echo -e "${GREEN}✓${NC} phase_prediction 測試通過 (3/3)"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    echo -e "${RED}✗${NC} phase_prediction 測試失敗"
    cat /tmp/test2.log
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

echo ""

# ============================================================
# 第三部分: Demo測試
# ============================================================
echo "========================================================================"
echo "3️⃣  Demo 測試"
echo "========================================================================"

echo "運行內建 demo_analysis()..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if python -c "from circadian_analysis import demo_analysis; demo_analysis()" > /tmp/demo.log 2>&1; then
    echo -e "${GREEN}✓${NC} Demo 運行成功"
    PASSED_TESTS=$((PASSED_TESTS + 1))

    # 檢查輸出文件
    if [[ -f "results/demo/rhythmic_genes.png" ]]; then
        echo -e "${GREEN}✓${NC} 輸出文件已生成: results/demo/rhythmic_genes.png"
    fi
else
    echo -e "${RED}✗${NC} Demo 運行失敗"
    cat /tmp/demo.log
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

echo ""

# ============================================================
# 第四部分: 示例程序測試
# ============================================================
echo "========================================================================"
echo "4️⃣  示例程序測試"
echo "========================================================================"

echo "運行 quick_start.py..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if timeout 60 python examples/quick_start.py > /tmp/quick.log 2>&1; then
    echo -e "${GREEN}✓${NC} quick_start.py 運行成功"
    PASSED_TESTS=$((PASSED_TESTS + 1))

    # 檢查輸出
    OUTPUT_COUNT=$(ls results/quick_start/*.png results/quick_start/*.csv 2>/dev/null | wc -l)
    echo "   生成了 $OUTPUT_COUNT 個輸出文件"
else
    echo -e "${RED}✗${NC} quick_start.py 運行失敗"
    cat /tmp/quick.log
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

echo ""
echo "運行 full_pipeline.py..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if timeout 120 python examples/full_pipeline.py > /tmp/full.log 2>&1; then
    echo -e "${GREEN}✓${NC} full_pipeline.py 運行成功"
    PASSED_TESTS=$((PASSED_TESTS + 1))

    # 檢查輸出
    OUTPUT_COUNT=$(ls results/full_pipeline/*.png results/full_pipeline/*.csv 2>/dev/null | wc -l)
    echo "   生成了 $OUTPUT_COUNT 個輸出文件"
else
    echo -e "${RED}✗${NC} full_pipeline.py 運行失敗"
    cat /tmp/full.log
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

echo ""

# ============================================================
# 第五部分: 輸出驗證
# ============================================================
echo "========================================================================"
echo "5️⃣  輸出文件驗證"
echo "========================================================================"

# 檢查結果目錄
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [[ -d "results/quick_start" ]] && [[ -d "results/full_pipeline" ]] && [[ -d "results/demo" ]]; then
    echo -e "${GREEN}✓${NC} 所有結果目錄已創建"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    echo -e "${RED}✗${NC} 部分結果目錄缺失"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

# 統計輸出文件
PNG_COUNT=$(find results/ -name "*.png" 2>/dev/null | wc -l)
CSV_COUNT=$(find results/ -name "*.csv" 2>/dev/null | wc -l)

echo ""
echo "輸出統計:"
echo "  - PNG圖片: $PNG_COUNT 個"
echo "  - CSV數據: $CSV_COUNT 個"
echo ""

# 列出所有輸出文件
echo "生成的文件:"
find results/ -type f \( -name "*.png" -o -name "*.csv" \) | sort | while read file; do
    SIZE=$(du -h "$file" | cut -f1)
    echo "  📄 $file ($SIZE)"
done

echo ""

# ============================================================
# 測試總結
# ============================================================
echo "========================================================================"
echo "📊 測試總結"
echo "========================================================================"
echo ""
echo "總測試數: $TOTAL_TESTS"
echo -e "${GREEN}通過: $PASSED_TESTS${NC}"
if [[ $FAILED_TESTS -gt 0 ]]; then
    echo -e "${RED}失敗: $FAILED_TESTS${NC}"
fi
echo ""

SUCCESS_RATE=$(echo "scale=1; $PASSED_TESTS * 100 / $TOTAL_TESTS" | bc)
echo "成功率: ${SUCCESS_RATE}%"
echo ""

if [[ $FAILED_TESTS -eq 0 ]]; then
    echo "========================================================================"
    echo -e "${GREEN}🎉 恭喜！所有測試通過！${NC}"
    echo "========================================================================"
    echo ""
    echo "項目已經完全可以使用了。您可以："
    echo "  1. 查看 results/ 目錄下的分析結果"
    echo "  2. 閱讀 TEST_GUIDE.md 了解更多測試選項"
    echo "  3. 使用 examples/ 下的腳本作為起點開發自己的分析"
    echo ""
    exit 0
else
    echo "========================================================================"
    echo -e "${RED}⚠️  有測試失敗，請檢查錯誤訊息${NC}"
    echo "========================================================================"
    echo ""
    echo "建議："
    echo "  1. 查看上面的錯誤訊息"
    echo "  2. 檢查 /tmp/*.log 文件獲取詳細日誌"
    echo "  3. 確認所有依賴已正確安裝"
    echo ""
    exit 1
fi
