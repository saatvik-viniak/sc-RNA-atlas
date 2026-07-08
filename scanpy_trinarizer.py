import numpy as np
from scipy.stats import beta as beta_dist
from tqdm.auto import tqdm
 
def trinarize(k, n, f=0.2, confidence_threshold=0.95, confidence_low=0.05):
    #same prior as they used previously
    alpha_prior = 1.5
    beta_prior = 1.5
    
    ## Now this is the bayesian update and here the posterior is evaluted
    alpha_post = k + alpha_prior
    beta_post = (n - k) + beta_prior
 
    #probability the true fraction expressing the gene is > f (20%) : 1 - cdf gives the area above f
    p = 1.0 - beta_dist.cdf(f, alpha_post, beta_post)
 
    #95% cutoff = on, 5% cutoff = off, anything in between is ambiguous
    calls = np.where(p >= confidence_threshold, 'on',
             np.where(p <= confidence_low, 'off', 'ambiguous'))
    return calls
 
def build_trinarization_matrix(adata, cluster_key='leiden', count_layer='raw_counts', gene_key='Gene'):
    counts = adata.layers[count_layer]    #full raw counts, cells x genes
    genes = adata.var[gene_key].astype(str).values   #gene symbols (var_names are ensembl ids)
    clusters = adata.obs[cluster_key].values
 
    results = {}
    for cluster in tqdm(np.unique(clusters), desc='trinarizing'):
        mask = clusters == cluster
        n = mask.sum() #total cells in the cluster
        k = np.asarray((counts[mask] > 0).sum(axis=0)).ravel() #cells with at least 1 count, per gene
        results[str(cluster)] = dict(zip(genes, trinarize(k, n)))
    return results
 