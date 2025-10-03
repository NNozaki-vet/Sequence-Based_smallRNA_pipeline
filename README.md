# Sequence-Based_smallRNA_pipeline
This repository hosts the original code used in the manuscript
“Sequence-based analysis reveals small-RNA profiles and their importance in transcriptome regulation in hepatocellular carcinoma across species” (paper link forthcoming).

The pipeline is for small RNA-seq (miRNA-seq). It retains exact read sequences, reduces reference-annotation bias, and captures small RNAs (e.g., tRNA, rRNA, YRNA, snoRNA, snRNA, and lncRNA fragments) in addition to canonical microRNA.

What’s inside in "Original_code"

1. Sequence filtering code

Filters and deduplicates reads while preserving exact sequences and counts.

This repository includes the code required for the filtering process used in our sequence-based analysis, scripts for generating evaluation data for the filtering procedure, and utilities for obtaining raw counts of the filtered sequences as well as for performing normalization (CPM and quantile normalization).


2. CLC Genomics Workbench annotation aggregation code

Utilities that import and aggregate CLC export tables into analysis-ready annotation matrices.

Before proceeding, run the CLC Genomics Workbench v22 (QIAGEN, Aarhus, Denmark) tools Quantify miRNA and Extract IsomiR Counts in sequence, and export the annotation data as CSV files into the specified folder.

More modules will be added if needed for manuscript reproduction (quantification, minimal stats, and figure helpers).
