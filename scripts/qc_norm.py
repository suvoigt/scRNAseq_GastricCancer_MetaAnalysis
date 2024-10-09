import os
import argparse
import numpy as np
import pandas as pd


import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix, issparse

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.vectors import IntVector, FloatVector, ListVector
import anndata


# activate automatic dataframe/array python to R conversion 
pandas2ri.activate()
numpy2ri.activate()

# import R packages
celda = importr('celda')
scDblFinder = importr('scDblFinder')
SingleCellExperiment = importr('SingleCellExperiment')
S4Vectors = importr('S4Vectors')
scuttle = importr('scuttle')
dropletutils = importr('DropletUtils')

# author: Susanne Voigt


# functions 

def _check_path(path):
    """
    Checks if the specified path exists as a file or directory.

    Args:
      path (str): The path to the file or directory.

    Returns:
      str: The validated path if it exists.

    Raises:
      argparse.ArgumentTypeError: If the file or directory path does not exist.
    """
    
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f"{path} does not exist.")

    return path


def load_data_csv(data, metadata, transpose=False, conversion=None):
    """
    Loads single-cell RNA-seq data and metadata, converting them into an AnnData object.

    Args:
        data (str): Path to the CSV file containing the raw count data.
        metadata (pd.DataFrame): DataFrame containing metadata with sample information.
        transpose (bool): Whether to transpose the data.
        conversion (str): Path to the conversion table file (optional).

    Returns:
        AnnData: The AnnData object containing the single-cell RNA-seq data with metadata.
    """

    if not os.path.exists(data):
        raise FileNotFoundError(f"Matrix file '{data}' does not exist.")
    
    # read data and convert to adata format
    if transpose:
        df = pd.read_csv(data, index_col=0).T
    else:
        df = pd.read_csv(data, index_col=0)

    adata = sc.AnnData(df)
    adata.obs_names = df.index
    adata.var_names = df.columns
    adata.var_names_make_unique()
    
    if conversion:
        conversion_table = pd.read_csv(conversion)
        mapping_dict = pd.Series(conversion_table.external_gene_name_grch38.values, index=conversion_table.external_gene_name_grch37).to_dict()
        converted_var_names = pd.Series(adata.var_names, index=adata.var_names).map(mapping_dict)
        adata = adata[:, converted_var_names.notna()]
        adata.var_names = converted_var_names.dropna().values
        adata.var_names_make_unique()

    adata.X = csr_matrix(adata.X)

    # add metadata
    meta = pd.DataFrame([metadata.loc[data].values] * adata.n_obs, columns=metadata.columns)
    meta.index = adata.obs_names
    adata.obs = meta

    return adata


def load_data_mtx(matrix, barcodes, features, metadata, transpose=False, filter_empty_drops=False, conversion=None):
    """
    Loads single-cell RNA-seq data from TSV and MTX files and metadata, converting them into an AnnData object.

    Args:
        matrix (str): Path to the matrix MTX file.
        barcodes (str): Path to the barcodes TSV file.
        features (str): Path to the features TSV file.
        metadata (pd.DataFrame): DataFrame containing metadata with sample information.
        transpose (bool): Whether to transpose the matrix.
        filter_empty_drops (bool): Whether to filter empty drops.
        conversion (str): Path to the conversion table file (optional).

    Returns:
        AnnData: The AnnData object containing the single-cell RNA-seq data with metadata.
    """

    if not os.path.exists(barcodes):
        raise FileNotFoundError(f"Barcodes file '{barcodes}' does not exist.")
    if not os.path.exists(features):
        raise FileNotFoundError(f"Features file '{features}' does not exist.")
    if not os.path.exists(matrix):
        raise FileNotFoundError(f"Matrix file '{matrix}' does not exist.")
    
    # read data 
    bar = pd.read_csv(barcodes, header=None, sep='\t').iloc[:, 0].values
    feat = pd.read_csv(features, header=None, sep='\t').iloc[:, 1].values

    if filter_empty_drops:

        if transpose:
            mat_raw = mmread(matrix).tocsr()
        else:
            mat_raw = mmread(matrix).T.tocsr()

        # prepare data for emptyDrops
        rows, cols = mat_raw.nonzero()
        data = mat_raw.data
        mat_r = ro.r['sparseMatrix'](i=IntVector(rows + 1), 
                                     j=IntVector(cols + 1),
                                     x=FloatVector(data),
                                     dims=IntVector(mat_raw.shape))
        # run emptyDrops
        ro.r('set.seed(42)')
        res = dropletutils.emptyDrops(mat_r)
        res_df = pandas2ri.rpy2py(ro.r['as.data.frame'](res))
        # filter out empty droplets based on FDR threshold
        cell_indices = np.where((res_df['FDR'] <= 0.01) & (~res_df['FDR'].isna()))[0]
        mat_filtered = mat_raw[:, cell_indices]        
        mat = mat_filtered.T.tocsr()
        bar = bar[cell_indices]

    else:

        if transpose:
            mat = mmread(matrix).T.tocsr()
        else:
            mat = mmread(matrix).tocsr()

    # get AnnData object
    adata = sc.AnnData(mat)
    adata.obs_names = bar
    adata.var_names = feat
    adata.var_names_make_unique()

    if conversion:
        conversion_table = pd.read_csv(conversion)
        mapping_dict = pd.Series(conversion_table.external_gene_name_grch38.values, index=conversion_table.external_gene_name_grch37).to_dict()
        converted_var_names = pd.Series(adata.var_names, index=adata.var_names).map(mapping_dict)        
        adata = adata[:, converted_var_names.notna()]
        adata.var_names = converted_var_names.dropna().values
        adata.var_names_make_unique()
    
    # add metadata
    meta = pd.DataFrame([metadata.loc[matrix].values] * adata.n_obs, columns=metadata.columns)
    meta.index = adata.obs_names
    adata.obs = meta

    return adata


def calculate_qc_metric(adata):
    """
    Calculates quality control (QC) metrics for an AnnData object.

    This function identifies mitochondrial, ribosomal, and hemoglobin genes in the dataset
    and calculates the proportions of counts attributed to these genes. It then computes
    various QC metrics and adds them to the AnnData object.

    Args:
        adata (AnnData): The AnnData object containing single-cell RNA-seq data.

    Returns:
        AnnData: The updated AnnData object with additional QC metrics.

    Note:
        This function relies on the `scanpy` library for QC metric calculation.
    """

    # calulate proportions of counts for mitochondrial, ribosomal and hemoglobin genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.str.contains(('^HB[^(P)]'))

    # calculate the QC metrics 
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo', 'hb'], 
                               inplace=True, percent_top=[20], log1p=True
                               )
    
    return adata


def filter_low_quality_cells(adata, min_cells=3, min_genes=500, max_genes=6000, min_counts=200, max_mito_pct=20):
    """
    Filters low-quality cells from an AnnData object based on several criteria.

    This function applies multiple filtering steps to remove low-quality cells and genes from the dataset.
    The filtering criteria include the minimum number of cells a gene must be expressed in, 
    the minimum and maximum number of genes expressed in a cell, the minimum total counts per cell,
    and the maximum percentage of mitochondrial genes.

    Args:
        adata (AnnData): The AnnData object containing single-cell RNA-seq data.
        min_cells (int): Minimum number of cells a gene must be expressed in to be retained. Default is 3.
        min_genes (int): Minimum number of genes a cell must express to be retained. Default is 500.
        max_genes (int): Maximum number of genes a cell can express to be retained. Default is 6000.
        min_counts (int): Minimum total counts per cell to be retained. Default is 200.
        max_mito_pct (float): Maximum percentage of mitochondrial gene counts for a cell to be retained. Default is 20.

    Returns:
        AnnData: The filtered AnnData object.
    """
    # filter data
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata = adata[adata.obs.pct_counts_mt < max_mito_pct, :].copy()

    return adata


def remove_ambient_rna(adata): 
    """
    Removes ambient RNA contamination from an AnnData object using the DecontX algorithm.

    This function prepares the input data for DecontX, runs the DecontX algorithm to 
    decontaminate the count matrix, and updates the AnnData object with both the original 
    and decontaminated counts stored in different layers.

    Args:
        adata (AnnData): The AnnData object containing single-cell RNA-seq data.

    Returns:
        AnnData: The updated AnnData object with ambient RNA removed.
    
     Note:
        This function relies on R and the `celda` package for DecontX, as well as the `rpy2` library 
        for interfacing between Python and R.
    """

    if issparse(adata.X):
        adata.X = adata.X.toarray()

    # prepare decontX input
    counts_r = numpy2ri.py2rpy(adata.X.T)
    # run decontX
    decontX_res = celda.decontX(counts_r)
    # extract decontaminated counts from the RS4 object
    decontX_counts_r = decontX_res.rx2('decontXcounts')
    # extract the data matrix from decontX_counts_r, transpose and round to integers
    decontX_counts_np = np.ceil(np.array(ro.r['as.matrix'](decontX_counts_r)).T)

    # store original and decontX counts as layers in AnnData object
    adata.layers['counts'] = csr_matrix(adata.X)
    adata.layers['decontX_counts'] = csr_matrix(decontX_counts_np)
    # overwrite original counts with decontX counts in .X
    adata.X = decontX_counts_np

    adata.X = csr_matrix(adata.X)

    return adata 


def detect_doublets(adata):
    """
    Detects doublets in an AnnData object using the scDblFinder method.

    This function prepares the input data for scDblFinder, runs the algorithm to detect doublets,
    and updates the AnnData object with doublet scores and classifications.

    Args:
        adata (AnnData): The AnnData object containing single-cell RNA-seq data.

    Returns:
        AnnData: The updated AnnData object with doublet scores and classifications.

    Note:
        This function relies on R and the `scDblFinder` package, as well as the `rpy2` library
        for interfacing between Python and R.
    """

    if issparse(adata.X):
        adata.X = adata.X.toarray()

    # prepare scDblFinder input 
    counts_r = numpy2ri.py2rpy(adata.X.T)
    sce_input = SingleCellExperiment.SingleCellExperiment(ListVector({'counts': counts_r}))
    # run scDblFinder
    ro.r('set.seed(42)')
    sce_results = scDblFinder.scDblFinder(sce_input)

    # extract results and add to adata 
    colData = sce_results.slots['colData']
    colData_df = pandas2ri.rpy2py(ro.r['as.data.frame'](colData))
    adata.obs['scDblFinder_score'] = colData_df['scDblFinder.score'].values
    adata.obs['scDblFinder_class'] = pd.Categorical(colData_df['scDblFinder.class'].values)

    adata.X = csr_matrix(adata.X)

    return adata


def shifted_log_norm(adata):
    """
    Applies shifted logarithmic normalization to the AnnData object.

    This function performs total-count normalization and then applies a log1p (logarithm of one plus the input) transformation.
    The results are stored in a new layer of the AnnData object. The function uses the `sc.pp.normalize_total` method 
    for total-count normalization. The `target_sum=None` parameter means  that each cell's counts are scaled to 
    the median count of the original counts per cell, preserving the total counts per cell rather than scaling to
    a fixed target sum.

    Args:
        adata (AnnData): The AnnData object containing single-cell RNA-seq data.

    Returns:
        AnnData: The AnnData object with an additional layer 'log1p_norm' containing the normalized and log-transformed data.
    """
    
    # normalize
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers['log1p_norm'] = sc.pp.log1p(scales_counts['X'], copy=True)


    return adata


def scran_norm(adata):
    """
    Applies scran normalization to the AnnData object.

    This function performs preliminary clustering and calculates size factors using the scran method,
    then normalizes the count data based on these size factors and stores the normalized data in a new layer of the AnnData object.

    Args:
        adata (AnnData): The AnnData object containing single-cell RNA-seq data.

    Returns:
        AnnData: The AnnData object with an additional layer 'scran_norm' containing the normalized and log-transformed data.
    """

    if issparse(adata.X):
        adata.X = adata.X.toarray()

    # preliminary clustering 
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added='groups', flavor='igraph', n_iterations=2)

    # get single cell object from count matrix
    sce = SingleCellExperiment.SingleCellExperiment(ListVector({'counts': adata_pp.X.T}))
    # calculate and size factors
    size_factors = scuttle.pooledSizeFactors(sce, clusters=adata_pp.obs['groups'])
    adata.obs['size_factors'] = size_factors
    # normalize 
    scran = adata.X / adata.obs['size_factors'].values[:, None]
    adata.layers['scran_norm'] = csr_matrix(sc.pp.log1p(scran))
    adata.X = csr_matrix(adata.X)

    return adata


if __name__ == "__main__":

    # define and parse command-line options:
    parser = argparse.ArgumentParser(description='Process single-cell RNA-seq data.')

    parser.add_argument('--data', type=_check_path, required=True, help='Path to raw data.')
    parser.add_argument('--metadata', type=_check_path, required=True, help='Path to .csv with metadata with first column including filenames.')
    parser.add_argument('--transpose', action='store_true', help='Flag to set transpose to True, if input matrix is of format genes x cells and needs to be transposed.')
    parser.add_argument('--filter_empty_drops', action='store_true', help='Flag to set filter_empty_drops to True, if input matrix needs to be filtered to remove empty droplets.')
    parser.add_argument('--conversion', type=_check_path, help='Path to the conversion table file (optional).')

    args = parser.parse_args()




    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--data', type=_check_path, required=True, help='Path to raw data.')
    parser.add_argument('--metadata', type=_check_path, required=True, help='Path to .csv with metadata with first column including filenames.')
    parser.add_argument('--transpose', action='store_true', help='Flag to set transpose to True, if input matrix is of format genes x cells and needs to be transposed.')
    parser.add_argument('--filter_empty_drops', action='store_true', help='Flag to set filter_empty_drops to True, if input matrix needs to be filtered to remove empty droplets.')
    parser.add_argument('--conversion', type=_check_path, help='Path to the conversion table file (optional).')
       

    args = parser.parse_args()

    # main
    os.chdir(args.data)
    metadata = pd.read_csv(args.metadata, index_col=0)

    for file in metadata.index.tolist():

        # load data
        print('Load raw data ...')
        if '.mtx' in file:
            adata = load_data_mtx(file, metadata.loc[file, 'barcodes_file'], metadata.loc[file, 'features_file'], metadata, args.transpose, args.filter_empty_drops, args.conversion)

        else:
            adata = load_data_csv(file, metadata, args.transpose, args.conversion) 

        # quality control
        print('Quality control ...')
        adata = calculate_qc_metric(adata)
        adata = filter_low_quality_cells(adata)
        adata = remove_ambient_rna(adata)
        adata = detect_doublets(adata)

        # normalize
        print('Normalize ...')
        adata= shifted_log_norm(adata)
        adata= scran_norm(adata)

        # save
        print('Save processed data ...')
        adata.write(f"{metadata.loc[file, 'dataset']}_{metadata.loc[file, 'sample_name']}_processed.h5ad", compression='gzip')


