#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# paramètres snake
tsv_file   <- snakemake@input[[1]]
pdf_file   <- snakemake@output[["pdf"]]
png_file   <- snakemake@output[["png"]]
padj_cut   <- snakemake@params[["padj"]]
lfc_cut    <- snakemake@params[["lfc"]]

# déduire les noms des deux conditions
contrast <- basename(tsv_file) |>
  sub("^deseq2_", "", x = _) |>
  sub("\\.tsv$",  "", x = _)
parts <- strsplit(contrast, "_vs_")[[1]]
cond1 <- parts[1];  cond2 <- parts[2]

# lecture des résultats DESeq2
dt <- read_tsv(tsv_file, show_col_types = FALSE) |>
  mutate(
    status = case_when(
      padj < padj_cut &  log2FoldChange >  lfc_cut ~ "Up",
      padj < padj_cut &  log2FoldChange < -lfc_cut ~ "Down",
      TRUE ~ "No"
    ),
    negLog10FDR = -log10(padj)
  )

# volcano plot
p <- ggplot(dt, aes(log2FoldChange, negLog10FDR, colour = status)) +
  geom_point(size = 1.2, alpha = .8) +
  scale_colour_manual(values = c(Down = "forestgreen",
                                 Up   = "firebrick",
                                 No   = "grey50")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cut),   linetype = "dashed") +
  labs(title = paste(cond1, "vs", cond2),
       x = expression(log[2]~Fold~Change),
       y = expression(-log[10]~FDR),
       colour = NULL) +
  theme_bw(base_size = 12)

# annoter 10 gènes plus significatifs
top10 <- dt |>
  filter(status != "No") |>
  arrange(padj) |>
  slice_head(n = 10)

p <- p + geom_text_repel(
  data = top10,
  aes(label = gene_id),
  size = 3,
  max.overlaps = 50
)
# calcul nombre total de gènes et annotation
total_feat <- nrow(dt)
caption    <- paste("Total:", total_feat, "features")

p <- p +
  annotate("text",
           x = max(dt$log2FoldChange, na.rm = TRUE) * 0.98,  
           y = min(dt$negLog10FDR,   na.rm = TRUE) + 0.05,   
           label = caption,
           hjust = 1, vjust = 0,
           size = 3)

# sauvegarde
dir.create(dirname(pdf_file), showWarnings = FALSE, recursive = TRUE)
ggsave(pdf_file, p, width = 7, height = 7)
ggsave(png_file, p, width = 7, height = 7, dpi = 300)
