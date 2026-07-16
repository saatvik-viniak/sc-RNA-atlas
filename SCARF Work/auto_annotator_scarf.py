# annotator.py
import glob
import os
import yaml


def load_rules(path):
    rules = []
    search_pattern = os.path.join(path, "*.md")

    for filename in glob.glob(search_pattern):
        with open(filename) as file:
            doc = next(yaml.load_all(file, Loader=yaml.SafeLoader))
        if doc is None:
            continue
        if "name" in doc:
            name = doc["name"]
        else:
            raise ValueError(
                os.path.basename(filename)
                + "this file didn't have a 'name' attribute which is required"
            )
        if "abbreviation" in doc:
            abbreviation = doc["abbreviation"]
        else:
            raise ValueError(
                os.path.basename(filename)
                + "this file didn't have a 'abbreviation' attribute which is required"
            )
        if "definition" in doc:
            definition = doc["definition"]
        else:
            raise ValueError(
                os.path.basename(filename)
                + "did not contain a 'definition' attribute, which is required."
            )
        if "categories" in doc:
            categories = doc["categories"]
        else:
            raise ValueError(
                os.path.basename(filename)
                + " did not contain a 'categories' attribute, which is required."
            )
        rules.append(
            {
                "name": name,
                "abbreviation": abbreviation,
                "definition": definition,
                "categories": categories,
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
        if all(g in present_genes for g in required) and all(
            g in absent_genes for g in excluded
        ):
            matches.append(rule["abbreviation"])
    return matches