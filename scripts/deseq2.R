#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)   # vroom
  library(dplyr)
})

## ── chemins fournis par Snakemake ──────────────────────────
counts_file <- snakemake@input[["counts"]]
meta_file   <- snakemake@input[["meta"]]
out_table   <- snakemake@output[["table"]]
out_rds     <- snakemake@output[["rds"]]

## ── lecture des données ────────────────────────────────────
cts  <- read_tsv(counts_file, col_types = cols())
meta <- read_tsv(meta_file,   col_types = cols()) |>
  mutate(across(where(is.character), factor))   # caractères → facteurs

# on place les gènes en rownames et on « arrondit » les comptes
mat <- as.matrix(round(cts[,-1]))
rownames(mat) <- cts$gene_id

## ── objet DESeq2 ───────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(mat, colData = meta, design = ~ condition)
dds <- DESeq(dds)

res <- results(dds) |>
  as_tibble(rownames = "gene_id") |>
  arrange(padj)

## ── sorties ────────────────────────────────────────────────
write_tsv(res,  out_table)
saveRDS(dds, out_rds)
