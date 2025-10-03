# Sequence-Based_smallRNA_pipeline
This repository hosts the original code used in the manuscript
“Sequence-based analysis reveals small-RNA profiles and their importance in transcriptome regulation in hepatocellular carcinoma across species” (paper link forthcoming).

The pipeline is for small RNA-seq (miRNA-seq). It retains exact read sequences, reduces reference-annotation bias, and captures small RNAs (e.g., tRNA, rRNA, YRNA, snoRNA, snRNA, and lncRNA fragments) in addition to canonical microRNA.

What’s inside

1.Sequence filtering code
Filters and deduplicates reads while preserving exact sequences and counts.

2.CLC Genomics Workbench annotation aggregation
Utilities that import and aggregate CLC export tables into analysis-ready annotation matrices.

More modules will be added if needed for manuscript reproduction (quantification, minimal stats, and figure helpers).
