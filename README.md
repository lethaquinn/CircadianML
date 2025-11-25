# 🧬 Circadian Transcriptomics Analysis Pipeline

**Advanced computational framework for discovering circadian rhythms and biomarkers in transcriptomic data**

## 🌟 Features
- 🔬 **JTK_CYCLE Algorithm**: Non-parametric circadian rhythm detection using Kendall’s tau correlation  
- 🤖 **ML Biomarker Discovery**: Machine-learning–based biomarker identification  
- 📊 **Advanced Visualisation**: Comprehensive plotting and dashboard capabilities  
- 🧬 **Real Data Integration**: Automated GEO database downloading and processing

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/lethaquinn/circadian-transcriptomics.git
cd circadian-transcriptomics

# Install dependencies
pip install -r requirements.txt

# Make the source package importable (for local development)
# Linux/macOS:
export PYTHONPATH=src
# Windows PowerShell:
# $env:PYTHONPATH = 'src'

# Run demo analysis (generates synthetic data and saves figures)
python -c "from circadian_analysis import demo_analysis; demo_analysis()"
```

The command above produces a synthetic circadian expression dataset, executes the full
rhythm detection and phase prediction stack, and stores summary figures in `results/demo/`.

Comprehensive usage examples are provided in
[`examples/quick_start.py`](examples/quick_start.py) and
[`examples/full_pipeline.py`](examples/full_pipeline.py).

## 🧱 Core Modules

| Area             | Module                                | Highlights                                                                 |
| ---------------- | ------------------------------------- | -------------------------------------------------------------------------- |
| Data Pipeline    | `circadian_analysis.data`             | GEO fetch helpers, batch-friendly normalisation                            |
| Rhythm Detection | `circadian_analysis.rhythm_detection` | JTK_CYCLE, Lomb–Scargle, wavelet spectra, cosinor regression               |
| Phase Prediction | `circadian_analysis.phase_prediction` | Random Forest, gradient boosting, MLP, ensemble predictors                 |
| Network Analysis | `circadian_analysis.network_analysis` | Co-expression graph construction, community discovery, phase coupling      |
| Visualisation    | `circadian_analysis.visualization`    | Phase wheel plots, ranked gene charts, network diagrams, dashboard bundler |

## 📊 Analysis Capabilities

* Circadian rhythm detection across multiple algorithms
* Biomarker discovery via classical ML and neural regressors
* Publication-ready plots and interactive dashboards
* Statistical validation via cross-validation metrics and harmonic summaries

## 🎯 Applications

* Circadian rhythm disruption analysis
* Clinical chronobiology research
* Drug chronotoxicity assessment
* Sleep deprivation biomarker discovery

## 🏆 Technical Highlights

* Mathematical rigour: Kendall’s tau correlation & cosinor regression
* Machine learning: Random Forests, gradient boosting, and MLP regressors
* Data integration: Support for GEO datasets and custom expression matrices
* Reproducibility: Deterministic algorithms with seed control

## 📚 Documentation

* Tutorial: Step-by-step analysis walkthrough
* API Reference: Complete function documentation
* Examples: Jupyter notebook demonstrations
* Methodology: Detailed algorithm explanations

## 📄 Citation

If you use this software in your research, please cite:

> Lethaquinn (2024). *Circadian Transcriptomics Analysis Pipeline: Advanced computational framework for circadian rhythm and biomarker discovery.* GitHub: [https://github.com/lethaquinn/circadian-transcriptomics](https://github.com/lethaquinn/circadian-transcriptomics)

## 📞 Contact

* Author: Lethaquinn
* Email: [2679066373@qq.com](mailto:2679066373@qq.com)
* GitHub: [@lethaquinn](https://github.com/lethaquinn)

> “Decoding the molecular rhythms of life through computational analysis.”
```








