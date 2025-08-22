# CircadianML
Detecting circadian rhythms in gene expression using ML | JTK_CYCLE, LSTM, Network Analysis | Application-ready pipeline
# CircadianML 🕐

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/lethaquinn/CircadianML?style=social)](https://github.com/lethaquinn/CircadianML)

## Machine Learning Pipeline for Circadian Transcriptome Analysis

A comprehensive Python pipeline for detecting, analyzing, and predicting circadian rhythms in gene expression data. This project bridges chronobiology and artificial intelligence to enable precision medicine applications.

### 🎯 Key Features

- **Rhythm Detection**: Implementation of JTK_CYCLE, Lomb-Scargle, and Cosinor analysis
- **ML Prediction**: Phase prediction using Random Forest, XGBoost, and LSTM models
- **Network Analysis**: Gene co-expression networks and circadian module detection
- **Clinical Applications**: Chronotype classification and jet lag simulation
- **Interactive Visualization**: Real-time dashboards and publication-ready figures

### 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/lethaquinn/CircadianML.git
cd CircadianML

# Install dependencies
pip install -r requirements.txt

# Run the pipeline
python src/main_pipeline.py --data GSE54650
