# README.md

This document provides a detailed description and usage instructions for the **single-cell RNA-seq data processing script**. This script is designed to process single-cell RNA-seq data, perform quality control, normalize the data, and detect doublets using various methods.

## Table of Contents
- [Purpose](#purpose)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [Functions](#functions)
  - [load_data_csv](#load_data_csv)
  - [load_data_mtx](#load_data_mtx)
  - [calculate_qc_metric](#calculate_qc_metric)
  - [filter_low_quality_cells](#filter_low_quality_cells)
  - [remove_ambient_rna](#remove_ambient_rna)
  - [detect_doublets](#detect_doublets)
  - [shifted_log_norm](#shifted_log_norm)
  - [scran_norm](#scran_norm)
- [Author](#author)

## Purpose

This script automates the **single-cell RNA sequencing (scRNA-seq)** data processing pipeline, including:
- Loading data (from CSV or MTX format).
- Quality control (QC) metrics computation.
- Filtering low-quality cells.
- Removing ambient RNA contamination using **DecontX**.
- Detecting doublets using **scDblFinder**.
- Normalizing the data using **shifted log normalization** and **scran normalization**.

The script integrates both Python and R libraries, allowing for seamless integration of the R-based methods like **DecontX** and **scDblFinder** with **AnnData** objects in Python.

## Usage

### Command-line Usage

```bash
python process_scRNA_data.py --data <path_to_data> --metadata <path_to_metadata.csv> [--transpose] [--filter_empty_drops] [--conversion <path_to_conversion_table.csv>]
```

- `--data`: Path to raw data file (either CSV or MTX).
- `--metadata`: Path to metadata CSV file containing sample information.
- `--transpose`: Flag indicating if the input matrix is transposed (genes x cells). Default is cells x genes.
- `--filter_empty_drops`: Flag to filter empty droplets using **DropletUtils::emptyDrops**.
- `--conversion`: Path to a conversion table CSV for gene name mapping (optional).

### Example

```bash
python process_scRNA_data.py --data /path/to/raw_data.mtx --metadata /path/to/metadata.csv --filter_empty_drops --conversion /path/to/conversion.csv
```

## Dependencies

- **Python 3.x**
  - `scanpy`
  - `pandas`
  - `numpy`
  - `scipy`
  - `rpy2` (for interfacing with R)
- **R**
  - `celda`
  - `scDblFinder`
  - `SingleCellExperiment`
  - `DropletUtils`
  - `scuttle`

Ensure the necessary Python and R libraries are installed before running the script.

## Functions

### `load_data_csv`
- **Purpose**: Loads single-cell RNA-seq data from CSV files and metadata, converting them into an AnnData object.
- **Arguments**:
  - `data`: Path to the CSV file containing raw count data.
  - `metadata`: DataFrame with metadata.
  - `transpose`: Whether to transpose the data matrix.
  - `conversion`: (Optional) Conversion table for gene names.
- **Returns**: An AnnData object.

### `load_data_mtx`
- **Purpose**: Loads scRNA-seq data from MTX (matrix), barcodes, and features files, and integrates metadata.
- **Arguments**:
  - `matrix`: Path to MTX file.
  - `barcodes`: Path to barcodes file.
  - `features`: Path to features file.
  - `metadata`: DataFrame with sample metadata.
  - `filter_empty_drops`: Flag to filter empty droplets.
- **Returns**: An AnnData object.

### `calculate_qc_metric`
- **Purpose**: Calculates QC metrics, including mitochondrial, ribosomal, and hemoglobin genes.
- **Arguments**:
  - `adata`: AnnData object containing the scRNA-seq data.
- **Returns**: Updated AnnData object with QC metrics.

### `filter_low_quality_cells`
- **Purpose**: Filters low-quality cells based on multiple criteria.
- **Arguments**:
  - `adata`: AnnData object.
  - `min_cells`: Minimum number of cells per gene (default = 3).
  - `min_genes`: Minimum genes per cell (default = 500).
  - `max_genes`: Maximum genes per cell (default = 6000).
  - `min_counts`: Minimum counts per cell (default = 200).
  - `max_mito_pct`: Maximum mitochondrial percentage (default = 20%).
- **Returns**: Filtered AnnData object.

### `remove_ambient_rna`
- **Purpose**: Removes ambient RNA contamination using **DecontX**.
- **Arguments**:
  - `adata`: AnnData object.
- **Returns**: AnnData object with decontaminated counts.

### `detect_doublets`
- **Purpose**: Detects doublets in the dataset using **scDblFinder**.
- **Arguments**:
  - `adata`: AnnData object.
- **Returns**: AnnData object with doublet scores and classifications.

### `shifted_log_norm`
- **Purpose**: Applies shifted logarithmic normalization.
- **Arguments**:
  - `adata`: AnnData object.
- **Returns**: AnnData object with normalized data in `log1p_norm` layer.

### `scran_norm`
- **Purpose**: Applies **scran** normalization to the AnnData object.
- **Arguments**:
  - `adata`: AnnData object.
- **Returns**: AnnData object with normalized data in `scran_norm` layer.

## Author
- **Susanne Voigt**

This script was developed to streamline the analysis and processing of single-cell RNA-seq datasets, integrating both Python and R tools for comprehensive data handling.