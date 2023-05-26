### import libraries

# standard libraries
import numpy as np
import pandas as pd

# math libraries
# import scipy.sparse.issparse
# import scipy.stats.ranksums
# import scipy.interpolate
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess as  sm_lowess
import math #to round up

# single cell libraries
import scanpy as sc
sc.settings.verbosity = 0 

### functions
def filter_data(
    adata, 
    mito_perc=5, 
    min_genes=700, 
    no_doublet=True, 
    no_negative=True,
):
    '''
    Filters raw count matrix.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    mito_perc : `int` (default: 5)
        Maximum percentage of mitochondrial genes in each individual cell.
    min_genes : `int` (default: 700)
        Minimum number of unique genes expressed in each individual cell. 
    no_doublet : `bool` (default: True)
        Remove doublets (as indicated by "Doublet" in adata.obs["hashtags"])
    no_negative : `bool` (default: True)
        Remove cells with no hashtag (as indicated by "Negative" in adata.obs["hashtags"])
    '''
    
    # annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    
    # filter cells with high % of mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs['pct_counts_mt'] < mito_perc, :]

    # filter cells with low no. of unique genes expressed
    sc.pp.filter_cells(adata, min_genes=min_genes)

    # filter out doublets and negatives
    if no_doublet==True:
        adata = adata[adata.obs["hashtags"]!= "Doublet",:]
    if no_negative==True:
        adata = adata[adata.obs["hashtags"]!= "Negative",:]
    
    # remove unnecessary obs and vars in anndata object
    if "hashtags" in adata.obs.columns:
        adata.obs = adata.obs[["hashtags"]]
    else:
        adata.obs = adata.obs[[]]
    adata.var = adata.var[[]]
    
    return(adata)

def pearson_residuals(
    counts, 
    theta=100,
):
    '''
    Computes analytical residuals for NB model with a fixed theta, 
    clipping outlier residuals to sqrt(N) as proposed in 
    Lause et al. 2021 https://doi.org/10.1186/s13059-021-02451-7
    
    Parameters
    ----------
    counts: `matrix` 
        Matrix (dense) with cells in rows and genes in columns
    theta: `int` (default: 100)
        Gene-shared overdispersion parameter
    '''
    
    counts_sum0 = np.sum(counts, axis=0)
    counts_sum1 = np.sum(counts, axis=1)
    counts_sum  = np.sum(counts)

    # get residuals
    mu = counts_sum1 @ counts_sum0 / counts_sum
    z = (counts - mu) / np.sqrt(mu + (np.square(mu)/theta))

    # clip to sqrt(n)
    n = counts.shape[0]
    z[z >  np.sqrt(n)] =  np.sqrt(n)
    z[z < -np.sqrt(n)] = -np.sqrt(n)

    return z

def get_hvgs(
    adata, 
    no_of_hvgs=2000, 
    theta=100,
):
    '''
    Function to select the top x highly variable genes (HVGs) 
    from an anndata object. 
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    no_of_hvgs: `int` (default: 2000)
        Number of hig
    theta: `int` (default: 100)
        Gene-shared overdispersion parameter used in pearson_residuals 
    '''
    
    # get pearson residuals
    if scipy.sparse.issparse(adata.X):
        residuals = pearson_residuals(adata.X.todense(),theta)
    else:
        residuals = pearson_residuals(adata.X,theta)
    
    # get variance of residuals
    residuals_variance = np.var(residuals,axis=0) 
    variances = pd.DataFrame({"variances":pd.Series(np.array(residuals_variance).flatten()), 
                              "genes":pd.Series(np.array(adata.var_names))})
    
    # get top x genes with highest variance
    hvgs = variances.sort_values(by="variances", ascending=False)[0:no_of_hvgs]["genes"].values
    
    return hvgs

def find_degs(
    adata1, 
    adata2, 
    top_x=50, 
    return_values="degs",
):
    '''
    Finds differentially expressed genes using Wilcoxon rank sum test on the 
    distribution of each gene in adata1 vs adata2.
    
    Parameters
    ----------
    adata1
        Annotated data matrix with all cells from selected cluster.
    adata2
        Annotated data matrix with all cells from all but selected cluster.
    top_x: `int` (default: 50)
        Number of DEGs to return
    return_values: `str` (default:"degs")
        String determing whether to return the degs or their p-values.
    '''
    
    ### find the genes that the sets have in common
    genes = adata1.var_names.intersection(adata2.var_names) 

    p_values = []
    
    ### get count matrices
    counts1 = pd.DataFrame(adata1.X.todense(), columns=adata1.var_names)
    counts2 = pd.DataFrame(adata2.X.todense(), columns=adata2.var_names)

    ### calculate p-values per gene
    for gene in genes:
        p_value = scipy.stats.ranksums(counts1[gene], counts2[gene])[1] #rank sum test
        p_values.append(p_value)
    p_values = pd.Series(p_values, index=adata1.var_names) #add gene names
    
    ### select top x DEGs
    p_values = p_values.sort_values()[0:top_x]
    degs = p_values.index
    
    if return_values == "degs":
        return(degs)
    
    elif return_values == "p_values":
        return(p_values)
    
    else:
        print("False entry for returned")
        
def calculate_geneset_scores(
    adata, 
    geneset,
):
    '''
    Scores each cell in the adata dataset for the marker genes of each
    celltype in the set of marker genes.
    
    Parameters
    ----------
    adata
        Annotated data matrix
    geneset
        Dataframe with in each column marker genes assigned to one
        specific cell type. The column names should include the 
        corresponding cell type names.
    '''
    ### create empty dataframe
    df_geneset_scores = pd.DataFrame(index=adata.obs.index) 
    
    ### loop over each cell type in the geneset
    for j in range(geneset.shape[1]):
        geneset_name = geneset.columns[j]
        sc.tl.score_genes(adata, list(geneset.iloc[:, j]), score_name="geneset_score")
        df_geneset_scores[geneset_name] = list(adata.obs["geneset_score"])
        del adata.obs["geneset_score"]
    
    return(df_geneset_scores)        

def integrate_datasets(adata, basis, batch_key="time", n_comps=100):
    '''
    Function to integrate multiple subsets in an Adata object using
    Scanorama (Hie, Bryson & Berger, 2019). The resulting Scanorama
    space is reduced in dimensions using PCA.
    
    Parameters
    ----------
    adata
        Annotated data matrix
    basis: `matrix` 
        Matrix (dense) with cells in rows and genes or reduced 
        dimensions in columns
    batch_key: `str` (default: "time")
        The name of the column in adata.obs that differentiates among 
        experiments/batches.Cells from the same batch must be 
        contiguously stored in adata
    n_comps: `int` (default: 100)
        Number of dimensions to return for the integrated space
    '''
    
    # add basis to adata
    adata.obsm["X_basis"] = basis
    
    # scanorama batch correction
    sc.external.pp.scanorama_integrate(adata, batch_key, basis="X_basis")
    
    # reduce dimensions 
    adata.obsm["X_scanorama_reduced"] = sc.tl.pca(adata.obsm["X_scanorama"], n_comps=n_comps)
    
    # create dataframe of scanorama space
    scanorama_df = pd.DataFrame(adata.obsm["X_scanorama_reduced"], index=adata.obs_names)
    
    return(scanorama_df)
        
def label_transfer(adata, batch_key="time", basis='X_scanorama_reduced', label_key="clusters",
                   reference="control", query="3h", no_neighbours=10):
    '''
    Function to transfer labels from a reference to a query. 
    Query and reference should both be included in one
    Anndata object. 
    
    Parameters
    ----------
    adata
        Annotated data matrix
    batch_key: `str` (default: "time")
        The name of the column in adata.obs that differentiates 
        reference from query
    basis: `str` (default: "X_scanorama_reduced")
        The name of the matrix in adata.obsm that is used to 
        calculate distance between cells
    label_key: `str` (default: "clusters")
        The name of the column in adata.obs which contains 
        the labels that have to be transferred
    reference: `str` (default: "control")
        The name that seperates the reference from the query
        in the adata.obs column indicated using batch_key
    query
        The name that seperates the query from the 
        reference in the adata.obs column indicated using 
        batch_key
    no_neighbours: `int` (default: 10)
        Number of neighbours to use for data integration
    '''

    distances = scipy.spatial.distance.cdist(adata[adata.obs[batch_key]==reference].obsm[basis],
                                             adata[adata.obs[batch_key]==query].obsm[basis], 
                                             metric='euclidean')
    df_distances = pd.DataFrame(distances,
                                index=adata[adata.obs[batch_key]==reference].obs[label_key], 
                                columns=adata[adata.obs[batch_key]==query].obs_names)
    neighbours = df_distances.apply(lambda x: pd.Series(x.nsmallest(no_neighbours).index))
    transferred_labels = neighbours.value_counts().idxmax()
    transferred_labels = pd.Series(transferred_labels, dtype="category")
    transferred_labels.index = adata[adata.obs[batch_key]==query].obs_names  
    
    return transferred_labels
    
def get_delta_expression(adata, genes, time_key="time", 
                         label_key="clusters", cluster=None):
    '''
    Computes the expression change for each timestep for one or 
    multiple genes. The expression change is either calculated 
    for the whole dataset or just for one cluster. 
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    genes: `array` 
        Array with all genes for which we want to get the 
        expression change in each timestep.
    time_key: `str` (default: "time")
        The name of the column in adata.obs that differentiates 
        between timepoints.
    label_key: `str` (default: "clusters")
        The name of the column in adata.obs that differentiates 
        between clusters.
    cluster: `str` (default: None)
        Which cluster to calculate the expression change for. 
        If not entered, the function will calculate expression
        change for the whole dataset.
    '''
    
    # get timepoints
    timepoints = adata.obs[time_key].cat.categories.values
    
    # get expression for all response genes in each cell
    if cluster is None:
        expression = adata[:,genes].X.toarray()
        expression = pd.DataFrame(expression,
                                  columns=genes,
                                  index=adata.obs[time_key])
    else:
        expression = adata[adata.obs[label_key]==cluster][:,genes].X.toarray()
        expression = pd.DataFrame(expression,
                                  columns=genes,
                                  index=adata[adata.obs[label_key]==cluster].obs[time_key])

    # calculate mean expression in each timepoint (for all response genes)
    mean_expr = expression.groupby(level=0).mean()

    # get delta expression in each timestep (for all response genes)
    delta_expr = mean_expr.diff()[1:len(timepoints)]

    return delta_expr

def get_expression(adata, genes, time_key="time",
                   label_key="clusters", cluster=None):
    '''
    Gets the expression in each timepoint for one or 
    multiple genes. The expression can be found
    for the whole dataset or just for one cluster. 
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    genes: `array` 
        Array with all genes for which we want to get the 
        expression change in each timestep.
    time_key: `str` (default: "time")
        The name of the column in adata.obs that differentiates 
        between timepoints.
    label_key: `str` (default: "clusters")
        The name of the column in adata.obs that differentiates 
        between clusters.
    cluster: `str` (default: None)
        Which cluster to calculate the expression change for. 
        If not entered, the function will find expression
        for the whole dataset.
    '''
    
    # get expression for all response genes in each cell
    if cluster is None:
        expression = adata[:,genes].X.toarray()
        expression = pd.DataFrame(expression,
                                  columns=genes,
                                  index=adata.obs[time_key])
    else:
        expression = adata[adata.obs[label_key]==cluster][:,genes].X.toarray()
        expression = pd.DataFrame(expression,
                                  columns=genes,
                                  index=adata[adata.obs[label_key]==cluster].obs[time_key])

    # calculate mean expression in each timepoint (for all response genes)
    mean_expr = expression.groupby(level=0).mean()

    return mean_expr
    
    

def bin_smooth(x, y, xgrid, sample_size=0.5, window_size=50):
    
    # take samples
    samples = np.random.choice(len(x), int(len(x)*sample_size), replace=True)
    x_sample = x[samples]
    y_sample = y[samples]
    
    if window_size >= len(samples):
            print("Unable to smooth: sample size is not bigger than window size.")
    
    # sort samples
    x_s_sorted = np.sort(x_sample)
    y_s_sorted = y_sample[np.argsort(x_sample)]

    window_halfsize = window_size/2
    
    y_smooth = []
            
    for idx, yi in enumerate(y_s_sorted):
        if window_halfsize < idx < (len(y_s_sorted) - window_halfsize):
            y_smooth.append(np.mean(y_s_sorted[int(idx-window_halfsize-1):int(idx+window_halfsize-1)]))
        else:
            y_smooth.append(None)
    
    # apply found funcction to x values on xgrid
    ygrid = scipy.interpolate.interp1d(x_s_sorted, y_smooth, fill_value='extrapolate')(xgrid) 

    return(ygrid)

def loess_smooth(x, y, xgrid, sample_size=0.5):
    
    # take samples
    samples = np.random.choice(len(x), int(len(x)*sample_size), replace=True)
    #samples = np.random.choice(len(x), 200, replace=True)
    x_sample = x[samples]
    y_sample = y[samples]
    
    y_smooth = sm_lowess(y_sample, x_sample, frac=0.4, it=5, return_sorted = False)
    
    # apply found function to x values on xgrid
    ygrid = scipy.interpolate.interp1d(x_sample, y_smooth, fill_value='extrapolate')(xgrid) 

    return(ygrid)

def bootstrap_smoothing(x, y, method="bin", sampling_rounds=20, sample_size=0.5, window_size=50):
    """
    This function takes x and y data (such as gene expression or gene score in 
    time) and smooths the datapoints. Method for smoothing is 'bin' or 'loess'.
    For bin smoothing a window size can be specified. 
    """
    
    xgrid = np.linspace(x.min(),x.max(),num=50) #set x-es that are compared in the end
    K = sampling_rounds #set number of samples for bootstrap
    
    if method == "bin":
        ygrids = np.stack([bin_smooth(x, y, xgrid, sample_size=sample_size, window_size=window_size) for k in range(K)]).T
        
    elif method == "loess":
        ygrids = np.stack([loess_smooth(x, y, xgrid, sample_size=sample_size) for k in range(K)]).T
        
    mean = np.nanmean(ygrids, axis=1)
    stderr = np.nanstd(ygrids, axis=1, ddof=0)
    
    return(xgrid, ygrids, mean, stderr)