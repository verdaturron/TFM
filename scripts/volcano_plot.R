#!/usr/bin/env Rscript
# volcano_plot.R — génère un volcano-plot PDF + PNG

suppressPackageStartupMessages({
  library(data.table)  
  library(ggplot2)
  library(ggrepel)
})

# Paramètres 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: volcano_plot.R <DE_results.tsv> <padj_cutoff> <lfc_cutoff> <outfile_base>")
}
de_file      <- args[1]
padj_cutoff  <- as.numeric(args[2])  
lfc_cutoff   <- as.numeric(args[3])   
outfile_base <- args[4]               


dt <- fread(de_file)


dt[, status := fifelse(padj < padj_cutoff & log2FoldChange >  lfc_cutoff, "Up",
                       fifelse(padj < padj_cutoff & log2FoldChange < -lfc_cutoff, "Down", "No"))]

## Pour ggplot, -log10(FDR)
dt[, neglog10FDR := -log10(padj)]

# Plot
p <- ggplot(dt, aes(x = log2FoldChange, y = neglog10FDR, colour = status)) +
  geom_point(size = 1.2, alpha = .8) +
  scale_colour_manual(values = c(Down = "forestgreen", Up = "firebrick", No = "grey50")) +
  geom_vline(xintercept =  c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff),           linetype = "dashed") +
  labs(title = "Differential Expression",
       x = expression(log[2]~Fold~Change),
       y = expression(-log[10]~FDR),
       colour = "Differential Expression") +
  theme_bw(base_size = 12)

# annoter les points les plus extrêmes
to_label <- dt[status != "No"][order(padj)][1:10]         
p <- p + geom_text_repel(data = to_label,
                         aes(label = gene_id),          
                         max.overlaps = 50,
                         size = 3)

# Export
ggsave(paste0(outfile_base, ".pdf"), p, width = 7, height = 7)
ggsave(paste0(outfile_base, ".png"), p, width = 7, height = 7, dpi = 300)
