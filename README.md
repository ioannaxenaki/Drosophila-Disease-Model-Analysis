# Drosophila-Disease-Model-Analysis
This project analyzes disease model annotations from the FlyBase database to explore the genetic links between Drosophila melanogaster (fruit flies) and human diseases. By studying conserved genes (orthologs) between flies and humans, researchers can test gene mutations and interactions in flies to predict disease mechanisms relevant to humans.   
This project integrates FlyBase data with Disease Ontology (DO) and the structured MONDO disease ontology to perform statistical analyses on these genetic relationships.

# Project Overview
The analysis is divided into sequential tasks that process the datasets, map ontology networks, and calculate statistical associations:
+ Data Extraction: Parses the FlyBase ```disease_model_annotations_fb_2025_02.tsv``` dataset, focusing on fly genes (FBgn IDs), Human Disease Ontology IDs (DO ID), and how the gene affects the disease (DO qualifiers such as "Ameliorates", "Exacerbates", or "Model of").
+ MONDO Ontology Graph Parsing: Loads ```mondo.json``` and parses its nodes and edges (specifically is_a relationships) to build a navigable tree of disease hierarchies.
+ Category Mapping: Maps all identified DO IDs to their corresponding MONDO IDs, and traverses the MONDO tree to extract the broad "Level 3" disease categories for each gene.
+ Statistical Analysis (Fisher's Exact Test): Builds 2x2 contingency tables for all unique combinations of Mondo Categories and DO qualifiers, and computes the two-sided Fisher's Exact Test to determine the p-value and odds ratio of the relationships.
+ Fold Change Calculation (Bonus 1): Calculates the expected value of genes for each combination (assuming independent events) and computes the observed vs. expected Fold Change.
+ Multiple Test Correction (Bonus 2): Applies Benjamini-Hochberg False Discovery Rate (FDR) multiple test correction (alpha=0.05) to rigorously identify relationships that reject the null hypothesis.

# Setup and Data 
To run this analysis, you need the following two datasets:
1. FlyBase Annotations: Download [```disease_model_annotations_fb_2025_02.tsv.gz```](https://flybase.org/releases/FB2025_02/precomputed_files/human_disease/disease_model_annotations_fb_2025_02.tsv.gz) from the FlyBase releases page and extract it.
2. MONDO Ontology: Download the structured [```mondo.json```](https://purl.obolibrary.org/obo/mondo.json) file from the OBO library.

# Usage
Execute the main script from your terminal:
```
python "Project_2&3.py"
```
