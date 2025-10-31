# GOseq Enrichment Analysis Pipeline

This document describes how to run the GOseq pipeline for bias‑corrected GO over‑representation analysis on SNP‑hit genes.

## Requirements

* R >= 4.5.0
* R packages: goseq, GenomicRanges, GenomicFeatures, dplyr, tidyr, stringr, BiocManager

Install missing packages:

```r
install.packages(c("goseq","dplyr","tidyr","stringr"))
if (!requireNamespace("GenomicFeatures", quietly=TRUE)) BiocManager::install("GenomicFeatures")
if (!requireNamespace("GenomicRanges", quietly=TRUE))  BiocManager::install("GenomicRanges")
```

## Input files

Place these files in your working directory:

* **GeneLengths\_RAD\_universe.txt**: two-column table with headers `gene_id` and `length` (integer exon lengths in bp for each RAD-covered gene).
* **Ityp\_GOterms.txt**: two-column table with headers `SeqName` (gene ID) and `GO.IDs` (semicolon-delimited GO terms, prefixed by namespace `F:`, `P:`, or `C:`).
* **SNPs\_with\_geneID\_all.tsv**: tab-delimited mapping of SNPs to genes, with at least columns `CHROM`, `POS`, `snpID`, and `gene_id` (NA for intergenic).
* **enrichment\_script.R**: the R script that executes the GOseq pipeline.

## Run pipeline

```r
setwd("/path/to/working/directory")

# Read gene lengths
gene_lengths <- read.delim("GeneLengths_RAD_universe.txt", header=TRUE, sep="\t")

# Read GO terms
go_tbl <- read.delim("Ityp_GOterms.txt", header=TRUE, sep="\t")

# Read SNP→gene mapping
 snp2gene <- read.delim("SNPs_with_geneID_all.tsv", header=TRUE, sep="\t")

# Define universe and phenotype vector
all_genes <- gene_lengths$gene_id
de_vec    <- as.integer(all_genes %in% unique(snp2gene$gene_id))
names(de_vec) <- all_genes

# Bias vector (gene lengths)
bias_vec <- gene_lengths$length[match(all_genes, gene_lengths$gene_id)]

# Build gene2cat list
gene2cat <- list()
for(i in seq_len(nrow(go_tbl))){
  gid <- go_tbl$SeqName[i]
  gos <- go_tbl$GO.IDs[i]
  if(gos==""||is.na(gos)) next
  terms <- strsplit(gos, ";")[[1]]
  terms <- sub("^[FPC]:","",terms)
  terms <- trimws(terms)
  gene2cat[[gid]] <- terms
}

# Load goseq
library(goseq)

# Estimate PWF and run enrichment
pwf <- nullp(DEgenes=de_vec, bias.data=bias_vec)
GO.wall <- goseq(pwf, gene2cat=gene2cat, method="Wallenius")

# Export results
write.table(GO.wall, "GOseq_all_terms.tsv", sep="\t", quote=FALSE, row.names=FALSE)
sigGO <- subset(GO.wall, over_represented_FDR<0.05)
write.table(sigGO, "GOseq_significant_terms.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```
