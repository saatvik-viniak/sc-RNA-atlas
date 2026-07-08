from scipy.sparse import issparse
import numpy as np
from scipy.stats import beta as beta_dist

def trinarize_gene_in_cluster(counts, f=0.2, confidence_threshold=0.95, confidence_low=0.05):
    counts = np.asarray(counts).ravel()
    n = len(counts)
    k = np.sum(counts > 0)

    alpha_post = k + 1.5
    beta_post = (n - k) + 2.0
    p_confidence = 1.0 - beta_dist.cdf(f, alpha_post, beta_post)

    if p_confidence >= confidence_threshold:
        call = "on"
    elif p_confidence <= confidence_low:
        call = "off"
    else:
        call = "ambiguous"

    return call, {"p_confidence": p_confidence, "k": k, "n": n, "frac": k / n if n else np.nan}

def build_trinarization_matrix(adata, cluster_key="leiden", count_layer="raw_counts", gene_key="Gene"):
    if count_layer not in adata.layers:
        raise ValueError(f"{count_layer} not found in adata.layers")

    matrix = adata.layers[count_layer]
    clusters = adata.obs[cluster_key].astype(str).values

    if gene_key in adata.var.columns:
        genes = adata.var[gene_key].astype(str).values
    else:
        genes = adata.var_names.astype(str).values

    results = {}
    stats = {}

    for cluster in np.unique(clusters):
        mask = clusters == cluster
        cluster_counts = matrix[mask]

        results[cluster] = {}
        stats[cluster] = {}

        for i, gene in enumerate(genes):
            gene_slice = cluster_counts[:, i]
            if issparse(gene_slice):
                gene_counts = gene_slice.toarray().ravel()
            else:
                gene_counts = np.asarray(gene_slice).ravel()

            call, score = trinarize_gene_in_cluster(gene_counts)
            results[cluster][gene] = call
            stats[cluster][gene] = score

    return results, stats
