##### ANALYSIS CODE #####
##### WRITTEN BY MATT JOGODNIK #####

## USER DEFINED PARAMETERS ##

# you will need to change these to the full paths on your machine
samples_csv = "/.../SASP2023/data/samples.csv"
STAR_dir = "/.../SASP2023/..."
qc_out = "/.../SASP2023/data/qc.csv"
RSEM_dir = "/.../SASP2023/..."
GSEA_expr_data = "/.../SASP2023/data/GSEA_counts.txt"

## IMPORTS ##

library(tidyverse)
library(tximport)
library(DESeq2)
library(biomaRt)
library(ggrepel)

# read in list of samples
samples <- read_csv(samples_csv)

## FUNCTIONS ##

# Given a DESeqDataSet object and a comparison find the difference
# in expression between the two conditions and convert ENSG to gene
# symbols
run_comparison <- function(dds, comparison) {
  res <- dds %>% 
    results(contrast=comparison) %>% 
    as.data.frame() %>%
    rownames_to_column("gene_id") %>% 
    as_tibble() %>% 
    arrange(desc(log2FoldChange))
  
  gene_key <- getBM(filters= "ensembl_gene_id",
                    attributes= c("ensembl_gene_id","external_gene_name"),
                    values=res$gene_id,mart=mart)
  res <- res %>% 
    inner_join(gene_key, by=c("gene_id"="ensembl_gene_id")) %>% 
    dplyr::select(gene_id, external_gene_name,log2FoldChange,pvalue,padj) %>% 
    dplyr::rename(
      gene_name = external_gene_name
    ) %>% 
    filter(gene_name != "") %>% 
    mutate(
      padj = replace_na(padj, 1)
    )
  return (res)
}

# Produce a volcano plot of differentially expressed genes

geom_volcano <- function(data, label_threshold, upper_limit) {
  data <- data %>% 
    mutate(
      sig = case_when(
        padj >= 0.05 ~ "ns",
        log2FoldChange > 1 ~ "up",
        log2FoldChange < -1 ~ "down",
        TRUE ~ "ns"
      )
    )
  
  upreg <- data %>% 
    filter(padj < 0.05 & log2FoldChange > label_threshold)
  
  data %>% 
    na.omit() %>% 
    ggplot(mapping=aes(x = log2FoldChange, y = -log10(padj), color=sig)) +
    geom_point(size=1) +
    scale_color_manual(values=c("#00539B", "grey30", "#A51C30")) +
    theme_light() +
    geom_text_repel(upreg, color = "#A51C30", mapping=aes(label=gene_name), size=4, max.overlaps=100,nudge_x=1, segment.alpha=0.25)+
    theme_classic(base_size=20) +
    labs(
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    ) +
    theme(legend.position="none") +
    geom_vline(xintercept=c(-1,1), color="grey50", linetype="dashed")+
    geom_hline(yintercept=-log10(0.05), color="grey50", linetype="dashed") +
    scale_y_continuous(expand=c(0, 0), limits=c(0, upper_limit)) +
    scale_x_continuous(breaks=c(-40,-35,-30,-25,-20,-15,-10,-5,-1,0,1,5,10,15,20,25,30,35,40))
}

## QUALITY CONTROL ##

# intialize data frames to store qc data from STAR files
df <- tibble()
df2 <- tibble()

file_path <- STAR_dir
files <- list.files(path = file_path, pattern="*.Log.final.out", full.names=TRUE)

# loop through log files and store number of reads to df
for (file in files) {
  log <- read_file(file) %>%
    str_split(pattern = "\n", simplify=TRUE) %>%  as.vector()
  cell_path <- str_split(file, ".Log")[[1]][1]
  cell <- str_split(cell_path, "STAR/")[[1]][2]
  num_reads <- log[6] %>% 
    str_split(pattern = "\t", simplify=TRUE) %>% as.vector()
  num_reads <- num_reads[2]
  
  df <- df %>% 
    bind_rows(tibble_row(
      cell = cell,
      num_reads = num_reads,
    ))
}

file_path <- STAR_dir
files <- list.files(path = file_path, pattern="*.ReadsPerGene.out.tab", full.names=TRUE)

# loop through reads files and store qc metrics to df2
for (file in files) {
  log <- read_file(file) %>%
    str_split(pattern = "\n", simplify=TRUE) %>%  as.vector()
  cell_path <- str_split(file, ".ReadsPerGene")[[1]][1]
  cell <- str_split(cell_path, "STAR/")[[1]][2]
  num_reads_unmapped <- log[1] %>% 
    str_split(pattern = "\t", simplify=TRUE) %>% as.vector()
  num_reads_unmapped <- num_reads_unmapped[2]
  num_reads_multimapped <- log[2] %>% 
    str_split(pattern = "\t", simplify=TRUE) %>% as.vector()
  num_reads_multimapped <- num_reads_multimapped[2]
  num_reads_nofeature <- log[3] %>% 
    str_split(pattern = "\t", simplify=TRUE) %>% as.vector()
  num_reads_nofeature <- num_reads_nofeature[2]
  num_reads_multfeatures <- log[4] %>% 
    str_split(pattern = "\t", simplify=TRUE) %>% as.vector()
  num_reads_multfeatures <- num_reads_multfeatures[2]
  
  df2 <- df2 %>% 
    bind_rows(tibble_row(
      cell = cell,
      num_reads_unmapped = num_reads_unmapped,
      num_reads_multimapped = num_reads_multimapped,
      num_reads_nofeature = num_reads_nofeature,
      num_reads_multfeatures = num_reads_multfeatures
    ))
}

# merge df and df2, output to qc_out
df %>% 
  inner_join(df2) %>% 
  write_csv(qc_out)

## IMPORT EXPRESSION DATA ##

# use tximport to collect RSEM files in format that DESeq2 can read
dir <- RSEM_dir
files <- file.path(dir, paste0(samples$sample, ".genes.results"))
names(files) <- samples$sample
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# prepare mapping of ENSG to gene symbols
mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror="useast")
gene_key <- getBM(filters= "ensembl_gene_id",
                  attributes= c("ensembl_gene_id","external_gene_name"),
                  values=rownames(txi.rsem$counts),mart=mart)

## DIFFERENTIAL EXPRESSION ##

# prepare table of conditions in format that DESeq2 needs
conditionTable <- samples %>% 
  dplyr::rename(
    id = sample
  ) %>% 
  dplyr::select(id, condition)

# correct genes with length 0 to 1 to prevent DESeq2 error
txi.rsem$length[txi.rsem$length == 0] <- 1

# import data to DESeq2, 2 data sets: 1 for heatmaps and 1 for GSEA
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, colData=conditionTable, design = ~ condition)
ddsGSEA <- DESeqDataSetFromTximport(txi.rsem, colData=conditionTable, design = ~ condition)

# run DESeq2
ddsTxi <- DESeq(ddsTxi)

# use DESeq2 to normalize counts separately for GSEA analysis
ddsGSEA <- estimateSizeFactors(ddsGSEA)
normalized_counts <- counts(ddsGSEA, normalized=TRUE)

## MCF10A HT-DNA VS. CONTROL ##

# collect results for comparison HT-DNA vs. control
MCF10A_DE_results <- ddsTxi %>% 
  run_comparison(c("condition", "MCF10A_ht", "MCF10A_control"))

# create HT-DNA custom gene set and export
MCF10A_DE_results %>% 
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  dplyr::select(gene_name) %>% 
  write_tsv("MCF10A_HTDNA_geneset.tsv", col_names=FALSE)

# export all results
MCF10A_DE_results %>% 
  write_csv("MCF10A_HT_DE_results.csv")

# create volcano plot with results, save
MCF10A_DE_results %>% 
  mutate(
    padj = if_else(-log10(padj) > 10, 10^-10, padj)
  ) %>% 
  dplyr::filter(log2FoldChange > -10) %>% 
  geom_volcano(8,10.5) +
  scale_y_continuous(n.breaks=6, expand=c(0, 0), limits=c(0,10.5)) +
  labs(
    title = "MCF10A HT-DNA vs. Control"
  )
ggsave("MCF10A_volcano.pdf", height=8, width=12)

# create heatmap
dev.new()
pdf("heatmap_mcf10a.pdf")

MCF10A_genes <- MCF10A_DE_results %>% 
  filter(padj < 0.05) %>% 
  pull(gene_id)

MCF10A_data <- txi.rsem$abundance[MCF10A_genes, colnames(txi.rsem$abundance[]) %in% c(
  "MCF10A_1_control",
  "MCF10A_2_ht",
  "MCF10A_3_control",
  "MCF10A_4_ht",
  "MCF10A_5_control",
  "MCF10A_6_ht")]

ht_mcf10a <- MCF10A_data %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  Heatmap(
    column_title = "MCF10A HT-DNA vs. Control",
    column_title_side = "top",
    width = unit(8, "cm"),
    name = "z-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    use_raster = FALSE,
    show_row_names = FALSE,
    column_order = c(
      "MCF10A_2_ht",
      "MCF10A_4_ht",
      "MCF10A_6_ht",
      "MCF10A_1_control",
      "MCF10A_3_control",
      "MCF10A_5_control"
    ),
    column_labels = c(
      "S23",
      "S24",
      "S25",
      "S26",
      "S27",
      "S28"
    )
  )

draw(ht_mcf10a)
dev.off()

## EXPORT NORMALIZED COUNTS FOR GSEA ##

# discard samples used to create custom gene set
gsea_samples <- samples %>% 
  filter(str_starts(
    condition, "MCF10A_"
  )) %>% 
  filter(condition != "MCF10A_ht" & condition != "MCF10A_control") %>% 
  pull(sample)

# prepare matrix in format for GSEA, convert ENSG to gene symbols, export
matrix_for_gsea <- normalized_counts %>% 
  as.data.frame() %>% 
  dplyr::select(gsea_samples) %>% 
  rownames_to_column("gene_id") %>% 
  inner_join(gene_key, by=c("gene_id" = "ensembl_gene_id")) %>% 
  dplyr::rename(
    Name = external_gene_name
  ) %>% 
  dplyr::select(-gene_id) %>% 
  filter(Name != "") %>% 
  mutate(Description = NA) %>% 
  dplyr::relocate(Name, Description)
matrix_for_gsea %>% 
  write_tsv(GSEA_expr_data)

# now, data can be run in the GSEA Java implementation