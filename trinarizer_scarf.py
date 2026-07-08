## new trinarizer
import pandas as pd
from scipy.stats import beta as beta_dist
import numpy as np
def trinarize_gene_in_cluster(counts, f=0.2, confidence_threshold=0.95, confidence_low=0.05):
    n = len(counts) #tells us the total cells in the cluster
    k = np.sum(counts > 0) #making sure only cells with at least 1 count
   
    #defining the same prior as they did previously
    alpha_prior = 1.5
    beta_prior = 2.0
    
    #does the bayesian update and evalutes the posterior
    alpha_post = k + alpha_prior
    beta_post = (n - k) + beta_prior
    
    #this calculates the probability that the true fraction of cells expressing gene x :is > f (20%)
    #the beta.cdf part integrates the function f and gives us the AUC 
    #then 1 - cdf gives the area above the function f
    p_confidence = 1.0 - beta_dist.cdf(f, alpha_post, beta_post)
    
    # now just put the 95% cutoff, so that means 95% confident, or FDR of 0.05, if doesn't fit into either, denote as ambigious
    if p_confidence >= confidence_threshold:
        call = 'on'
    elif p_confidence <= confidence_low:
        call = 'off'
    else:
        call = 'ambiguous' 
        
    scores = {'p_confidence': p_confidence}
    return call, scores

def build_trinarization_matrix(ds, cluster_key='RNA_leiden_cluster'): # uses the leiden clusters for the # of clusters/their coordinates - can be switched to your choice
    x = ds.cells.fetch_all(cluster_key)
    genes = ds.RNA.feats.fetch_all('names')  # just fetches everything from scarf (all of the GENE names)
    clusters = np.unique(x) 
    # create the results dictionary here
    results = {}
    for cluster in clusters:
        mask = x == cluster 
        cluster_counts = ds.RNA.rawData[mask].compute() 
        results[cluster] = {}
        for i, gene in enumerate(genes):
            gene_counts = cluster_counts[:,i].flatten() 
            call, _ = trinarize_gene_in_cluster(gene_counts)
            results[cluster][gene] = call
    return results
