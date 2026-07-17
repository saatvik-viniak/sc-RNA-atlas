## scanpy auto annotator
import glob
import os
import yaml


def load_rules(path):
    rules = []
    for filename in glob.glob(os.path.join(path, "*.md")):
        with open(filename) as f:
            doc = next(yaml.load_all(f, Loader=yaml.SafeLoader))
        if doc is None:
            continue
        for key in ("name", "abbreviation", "definition", "categories"):
            if key not in doc:
                raise ValueError(f"{os.path.basename(filename)} is missing required attribute '{key}'")
        rules.append({k: doc[k] for k in ("name", "abbreviation", "definition", "categories")})
    return rules


def parse_definition(definition_str):
    genes = definition_str.strip().split()
    required = [x[1:] for x in genes if x.startswith("+")]
    excluded = [x[1:] for x in genes if x.startswith("-")]
    return required, excluded


def annotate_cluster(probs, rules):
    matches = []
    for rule in rules:
        required, excluded = parse_definition(rule["definition"])
        #start at 1 and multiply: p for each positive marker, (1 - p) for each negative marker
        p = 1.0
        for gene in required:
            p = p * probs.get(gene, 0.0)
        for gene in excluded:
            p = p * (1.0 - probs.get(gene, 0.0))
        if p > 0.5: #tag applies when the combined probability is greater than 50%
            matches.append(rule["abbreviation"])
    return matches


def auto_annotate_results(results, rules):
    ## results is {cluster: {gene: probability}} from build_trinarization_matrix
    return {cluster: annotate_cluster(probs, rules) for cluster, probs in results.items()}


def add_annotations_to_adata(adata, annotations, cluster_key="leiden", obs_key="auto_annotation"):
    ## maps annotation list back to each cell via its cluster label
    adata.obs[obs_key] = adata.obs[cluster_key].astype(str).map(
        lambda c: ", ".join(annotations.get(c, ["unassigned"]))
    )
    return adata