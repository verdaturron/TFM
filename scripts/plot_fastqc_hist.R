#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

raw_zip  <- snakemake@input[["raw_zip"]]
trim_zip <- snakemake@input[["trim_zip"]]
out_pdf  <- snakemake@output[["pdf"]]

# ────────────────── fonctions utilitaires ──────────────────
extract_module <- function(zip_path, module) {
  inner <- unzip(zip_path, list = TRUE)$Name
  inner <- inner[grepl("fastqc_data.txt$", inner)][1]
  txt   <- readLines(unz(zip_path, inner))
  start <- grep(paste0(">>", module), txt)
  end   <- grep(">>END_MODULE", txt); end <- end[end > start][1]
  block <- txt[(start + 1):(end - 1)]
  block <- block[!startsWith(block, "#")]
  read_tsv(paste(block, collapse = "\n"), show_col_types = FALSE)
}

read_len  <- function(z, lbl) extract_module(z, "Sequence Length Distribution") |>
  setNames(c("Length", "Count")) |>
  mutate(source = lbl)
read_qual <- function(z, lbl) extract_module(z, "Per sequence quality scores")  |>
  setNames(c("Quality", "Count")) |>
  mutate(source = lbl)

# ────────────────── données ─────────────────────────────────
len_all <- bind_rows(read_len(raw_zip,  "Before"),
                     read_len(trim_zip, "After")) |>
  mutate(Length = as.numeric(sub("^([0-9]+).*", "\\1", Length)) + 0.5)

qual_all <- bind_rows(read_qual(raw_zip,  "Before"),
                      read_qual(trim_zip, "After"))

# ────────────────── graphiques ──────────────────────────────
p_len <- ggplot(len_all, aes(Length, Count, fill = source)) +
  geom_col(alpha = .6, position = "identity") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  scale_fill_manual(values = c(Before = "grey60",
                               After  = "steelblue")) +
  theme_bw(11) +
  labs(title = "Read length", x = "bp", y = "Count")

p_qual <- ggplot(qual_all, aes(Quality, Count, colour = source)) +
  geom_step(linewidth = .7, direction = "mid") +
  scale_colour_manual(values = c(Before = "grey50",
                                 After  = "steelblue")) +
  theme_bw(11) +
  labs(title = "Per-sequence quality", x = "Phred", y = "Count")

# ────────────────── sortie ──────────────────────────────────
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf,
       gridExtra::arrangeGrob(p_len, p_qual, ncol = 2),
       width = 9, height = 4)
