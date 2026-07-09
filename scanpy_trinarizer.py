## scanpy trinarizer
import numpy as np
from scipy.stats import beta as beta_dist
from tqdm.auto import tqdm

def trinarize(k, n, f=0.2):
    #same prior as they used previously
    alpha_post = k + 1.5
    beta_post = (n - k) + 2.0

    #probability the true fraction expressing the gene is > f (20%); 1 - cdf gives the area above f
    return 1.0 - beta_dist.cdf(f, alpha_post, beta_post)

def build_trinarization_matrix(adata, cluster_key='leiden', count_layer='raw_counts', gene_key='Gene'):
    counts = adata.layers[count_layer]   #full raw counts, cells x genes
    genes = adata.var[gene_key].astype(str).values #gene symbols (var_names are ensembl ids)
    clusters = adata.obs[cluster_key].values

    results = {}
    for cluster in tqdm(np.unique(clusters), desc='trinarizing'):
        mask = np.asarray(clusters == cluster)              
        n = mask.sum() #total cells in the cluster
        k = np.asarray((counts[mask] > 0).sum(axis=0)).ravel()   #cells with at least 1 count, per gene
        results[str(cluster)] = dict(zip(genes, trinarize(k, n))) #store the continuous probability per gene as they do
    return results

