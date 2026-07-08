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


def annotate_cluster(present_genes, absent_genes, rules):
    matches = []
    for rule in rules:
        required, excluded = parse_definition(rule["definition"])
        if all(g in present_genes for g in required) and all(g in absent_genes for g in excluded):
            matches.append(rule["abbreviation"])
    return matches


def auto_annotate_results(results, rules):
    ## results is {cluster: {gene: 'on'/'off'/'ambiguous'}} from build_trinarization_matrix
    annotations = {}
    for cluster, calls in results.items():
        present_genes = {gene for gene, call in calls.items() if call == "on"}
        absent_genes  = {gene for gene, call in calls.items() if call == "off"}
        annotations[cluster] = annotate_cluster(present_genes, absent_genes, rules)
    return annotations


def add_annotations_to_adata(adata, annotations, cluster_key="leiden", obs_key="auto_annotation"):
    ## maps annotation list back to each cell via its cluster label
    adata.obs[obs_key] = adata.obs[cluster_key].astype(str).map(
        lambda c: ", ".join(annotations.get(c, ["unassigned"]))
    )
    return adata