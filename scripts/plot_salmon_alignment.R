#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite); library(dplyr); library(ggplot2)
})

out_pdf <- snakemake@output[["pdf"]]
jsons   <- snakemake@input          # liste des meta_info.json

rates <- lapply(jsons, function(j) {
  meta <- read_json(j)
  tibble(
    sample  = basename(dirname(dirname(j))),   # <- remonte d’un niveau de plus
    aligned = 100 * meta$num_mapped / meta$num_processed
  )
}) %>% bind_rows()

p <- ggplot(rates, aes(sample, aligned)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", aligned)), vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .05))) +
  labs(title = "Salmon – alignment rate per sample",
       x = NULL, y = "Aligned reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf, p, width = 8, height = 4.5)
message("✅ Barplot alignment écrit → ", out_pdf)
