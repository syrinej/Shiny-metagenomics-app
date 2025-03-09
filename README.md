# 📊 Metagenomics Data Visualization Shiny App

An interactive Shiny web application developed for visualizing and analyzing metagenomics data using the `phyloseq` R package.

##  Features

- Interactive visualizations:
  - Bar Plot (Taxonomic Composition)
  - Heatmap (Taxonomic Abundance)
  - Alpha Diversity Metrics (Shannon, Simpson)

- Easy-to-use interface for uploading custom datasets (OTU table, taxonomy, and sample metadata).


## Getting Started

### Prerequisites
- R (≥4.0 recommended)
- RStudio (recommended)

### Installation
Install required packages in R:
```R
install.packages("shiny")
install.packages("phyloseq")
install.packages("ggplot2")
