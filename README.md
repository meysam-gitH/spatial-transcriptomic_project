# Spatial Transcriptomics Analysis (Scanpy + Squidpy)

## Overview
This project implements an end-to-end pipeline for spatial transcriptomics data analysis using Scanpy and Squidpy.

## Features
- Quality control (QC)
- Normalization and log transformation
- PCA, UMAP, Leiden clustering
- Marker gene analysis
- Cell type annotation
- Spatial visualization

## Dataset
This project uses a sample Visium dataset provided by Scanpy for demonstration purposes.

## Project Structure
src/ # pipeline code
results/ # generated figures
data/ # input data (not included)

## How to Run
```bash
python main.py

Output
UMAP plot of clusters
Spatial distribution of cell types
