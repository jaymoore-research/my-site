# 🧠 Brain Atlas Viewer

Interactive single-cell RNA-seq data exploration dashboard.

## Data
- **Cells:** 59,505
- **Genes:** 58,232
- **Source:** [Allen Institute - SEA-AD: Seattle Alzheimer's Disease Brain Cell Atlas](https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443)

## Features

### 📊 Visualizations
- **UMAP/t-SNE plots** - Interactive embeddings colored by metadata or gene expression
- **Violin plots** - Gene expression across cell types
- **Cell statistics** - Counts, QC metrics, distributions

### 🎨 Coloring Options
- Cell type
- Disease status
- Tissue
- Development stage
- Gene expression (any gene)

### 🔍 Filters
- Cell type
- Disease
- Tissue

### 💾 Export
- Download metadata as CSV

## Installation

```bash
cd /Users/jaymoore/Desktop/brain-atlas
pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

The dashboard will open at `http://localhost:8501`

## Quick Start

1. Launch the app
2. Use the sidebar to select embeddings (UMAP/t-SNE)
3. Choose what to color by (cell type, disease, or gene expression)
4. Apply filters to focus on specific cell populations
5. Search for genes to see expression patterns
6. Explore violin plots for detailed expression analysis