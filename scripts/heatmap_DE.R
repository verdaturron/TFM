#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)       
  library(dplyr)       
  library(tidyr)       
  library(pheatmap)    
  library(tibble)})


de_file   <- snakemake@input[["de"]]      
counts_tsv<- snakemake@input[["counts"]]  
meta_tsv  <- snakemake@input[["meta"]]    
out_pdf   <- snakemake@output[["pdf"]]
out_png   <- snakemake@output[["png"]]
top_n     <- as.integer(snakemake@params[["top"]])  

# 1 gènes DE les plus significatifs
de <- read_tsv(de_file, show_col_types = FALSE) |>
  arrange(padj) |>
  slice_head(n = top_n) |>
  select(gene_id)

# 2 matrice d’expression normalisée 
cts <- read_tsv(counts_tsv, show_col_types = FALSE)
mat <- cts |> filter(gene_id %in% de$gene_id) |>
  column_to_rownames("gene_id") |>
  as.matrix()

# log-TPM 
mat <- log2(mat + 1)

# 3 annotation des échantillons
meta <- read_tsv(meta_tsv, show_col_types = FALSE) |>
  select(sample, condition) |>
  column_to_rownames("sample")

# 4 Heatmap
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

for (f in c(out_pdf, out_png)) {
  if (grepl("\\.png$", f)) png(f, 1200, 1400, res = 150)
  else                     pdf(f, 8, 10)
  
  pheatmap(mat,
           annotation_col = meta,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           fontsize_row = 6,
           main = basename(de_file))
  
  dev.off()
}

message("✓ Heatmap écrite → ", out_pdf, " + .png")
