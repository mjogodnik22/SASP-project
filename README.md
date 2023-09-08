# Project Repository

Analysis and code by Matt Jogodnik
bulkRNAseq pipeline adapted from Nikos Mynhier

This repository contains the code required to run the analysis contained in the article.

## src

### pipeline

The pipeline folder contains the snakemake pipeline used to map fastqs to the reference genome and quantify expression levels. Modify the paths in Config.yaml as needed. pipeline_environment.yml contains the specifications for the conda environment used.

### analysis.R

analysis.R contains all the code used to generate a normalized counts matrix and custom gene set as input for GSEA. Modify the file paths as needed.

## data

The data folder contains all the files needed to run GSEA analysis including:

- The custom gene set file (MCF10A_HTDNA_geneset.gmx)
- The sample condition mapping file (samples.cls)
- The normalized counts matrix - **generate first from analysis.R** (GSEA_counts.txt)

### MCF10A_HTDNA_geneset.gmx

This file contains the genes belonging to the custom gene set to be used in GSEA. This can be easily re-generated from the analysis.R output *MCF10A_HTDNA_geneset.tsv*.

### samples.cls

This file maps each sample to its experimental condition. See GSEA documentation for proper formatting. Ensure the mappings match the order of columns in *GSEA_counts.txt*.