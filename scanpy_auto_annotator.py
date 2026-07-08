import os
import yaml


def load_rules(path):
    rules = []

    for root, _, files in os.walk(path):
        for file in files:
            if file == "README.md":
                continue
            if not (file.endswith(".md") or file.endswith(".yaml")):
                continue

            filename = os.path.join(root, file)

            with open(filename) as f:
                doc = next(yaml.load_all(f, Loader=yaml.SafeLoader))

            if doc is None:
                continue

            if "name" not in doc:
                raise ValueError(os.path.basename(filename) + " did not contain a 'name' attribute, which is required.")
            if "abbreviation" not in doc:
                raise ValueError(os.path.basename(filename) + " did not contain an 'abbreviation' attribute, which is required.")
            if "definition" not in doc:
                raise ValueError(os.path.basename(filename) + " did not contain a 'definition' attribute, which is required.")

            rules.append(
                {
                    "name": doc["name"],
                    "abbreviation": doc["abbreviation"],
                    "definition": doc["definition"],
                    "categories": doc.get("categories", None),
                }
            )

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
    annotations = {}

    for cluster, gene_calls in results.items():
        present_genes = {gene for gene, call in gene_calls.items() if call == "on"}
        absent_genes = {gene for gene, call in gene_calls.items() if call == "off"}
        annotations[cluster] = annotate_cluster(present_genes, absent_genes, rules)

    return annotations


def add_annotations_to_adata(adata, annotations, cluster_key="leiden", obs_key="auto_annotation"):
    cluster_to_label = {
        cluster: " ".join(labels) if len(labels) > 0 else ""
        for cluster, labels in annotations.items()
    }
    adata.obs[obs_key] = adata.obs[cluster_key].map(cluster_to_label)
    return adata
