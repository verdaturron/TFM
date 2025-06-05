#!/usr/bin/env Rscript
library(tximport)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# arguments Snakemake

quant_files <- unlist(snakemake@input[["quant"]]) 
names(quant_files) <- basename(dirname(quant_files))


gtf_file   <- snakemake@input[["gtf"]]
out_counts <- snakemake@output[["counts"]]
out_tpm    <- snakemake@output[["tpm"]]

# table tx2gene 
tx2gene <- read_tsv(gtf_file,
                    comment = "#", col_names = FALSE,
                    col_types = cols(.default = "c")) |>
  filter(X3 == "transcript") |>
  mutate(
    tx_id   = str_match(X9, 'transcript_id "([^"]+)"')[, 2],
    gene_id = str_match(X9, 'gene_id "([^"]+)"')[, 2]
  ) |>
  select(tx_id, gene_id) |>
  drop_na() |>
  distinct()

# tximport
txi <- tximport(files   = quant_files,
                type    = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

write_tsv(as_tibble(txi$counts,    rownames = "gene_id"), out_counts)
write_tsv(as_tibble(txi$abundance, rownames = "gene_id"), out_tpm)
