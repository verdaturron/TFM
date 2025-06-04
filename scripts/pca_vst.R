#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)   # déjà présent dans l’env heatmap.yaml
})

# ----- I/O -----
cts_file  <- snakemake@input[["counts"]]
meta_file <- snakemake@input[["meta"]]
out_pdf   <- snakemake@output[["pdf"]]
out_png   <- snakemake@output[["png"]]

# ----- data -----
cts  <- readr::read_tsv(cts_file, show_col_types = FALSE)
meta <- readr::read_tsv(meta_file, show_col_types = FALSE) |>
  mutate(across(where(is.character), factor)) |>
  column_to_rownames("sample")

mat  <- as.matrix(round(cts[,-1]))
rownames(mat) <- cts$gene_id

# même ordre échantillons
mat  <- mat[, rownames(meta)]

# ----- VST + PCA -----
dds     <- DESeqDataSetFromMatrix(mat, meta, design = ~ condition)
vsd     <- vst(dds, blind = TRUE)
pca     <- prcomp(t(assay(vsd)))
percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pdat <- as.data.frame(pca$x[,1:2]) |>
  rownames_to_column("sample") |>
  left_join(meta |> rownames_to_column("sample"), by = "sample")

p <- ggplot(pdat, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3, alpha = .9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50, show.legend = FALSE) +
  labs(title = "PCA on VST counts",
       x = paste0("PC1 (", percent[1], "%)"),
       y = paste0("PC2 (", percent[2], "%)")) +
  theme_bw(12) +
  theme(panel.grid = element_blank(),
        legend.position = "right")

# ----- write -----
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf, p, width = 6, height = 5)
ggsave(out_png, p, width = 6, height = 5, dpi = 300)
message("✅ PCA plot written → ", out_pdf)
