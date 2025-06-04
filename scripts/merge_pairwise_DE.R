#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)    # read_tsv / write_tsv
  library(dplyr)    # bind_rows
  library(stringr)  # str_remove
})

## ── chemins fournis par Snakemake ──────────────────────────────────
outfile <- snakemake@output[[1]]   # results/deseq2/deseq2_all_pairwise.tsv
inputs  <- snakemake@input         # vecteur de 6 fichiers TSV pairwise

## ── lecture + ajout de la colonne 'contrast' ───────────────────────
all_tbl <- lapply(inputs, function(f) {
  contrast <- f %>% basename() %>%
    str_remove("^deseq2_") %>% str_remove("\\.tsv$")
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(contrast = contrast, .before = 1)
})

## ── fusion des tables ──────────────────────────────────────────────
merged_tbl <- bind_rows(all_tbl)

## ── écriture ───────────────────────────────────────────────────────
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
write_tsv(merged_tbl, outfile)
message("✅ Écrit : ", outfile)
