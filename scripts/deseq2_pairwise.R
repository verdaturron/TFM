#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# ────────────────────────────────────────────────────────────
# Chemins passés par Snakemake
# ────────────────────────────────────────────────────────────
counts_file <- snakemake@input[["counts"]]      # results/expression/gene_counts.tsv
meta_file   <- snakemake@input[["meta"]]        # config/samples.tsv
cond1       <- snakemake@wildcards[["cond1"]]
cond2       <- snakemake@wildcards[["cond2"]]
out_tsv     <- snakemake@output[["table"]]      # results/deseq2/deseq2_<cond1>_vs_<cond2>.tsv

# ────────────────────────────────────────────────────────────
# Lecture du fichier counts
# ────────────────────────────────────────────────────────────
cts_df <- read_tsv(counts_file, col_types = cols())

# Nettoyage du header : supprime les espaces invisibles
clean_names <- function(x) str_trim(x, side = "both")
colnames(cts_df) <- clean_names(colnames(cts_df))

# Noms des échantillons (sans la 1re colonne gene_id)
counts_names <- colnames(cts_df)[-1]

# ────────────────────────────────────────────────────────────
# Lecture du fichier metadata
# ────────────────────────────────────────────────────────────
meta <- read_tsv(meta_file, col_types = cols()) %>%
  mutate(across(where(is.character), factor)) %>%
  # nettoie les colonnes
  mutate(sample    = clean_names(sample),
         condition = clean_names(condition)) %>%
  column_to_rownames("sample")

meta_names <- rownames(meta)

# ── Vérification concordance jeu de counts ↔ metadata ───────
if (!setequal(counts_names, meta_names)) {
  stop(
    "❌ Les échantillons de gene_counts.tsv et de samples.tsv ne correspondent pas.\n",
    "   counts : [", paste(counts_names, collapse = ", "), "]\n",
    "   meta   : [", paste(meta_names,   collapse = ", "), "]\n"
  )
}

# Réordonne meta pour qu’il suive EXACTEMENT l’ordre des colonnes counts
meta <- meta[counts_names, , drop = FALSE]

# Double-vérification
if (!identical(rownames(meta), counts_names)) {
  stop("❌ Les noms d’échantillon sont les mêmes mais ne sont toujours pas dans le même ordre.")
}

# ────────────────────────────────────────────────────────────
# Matrice de comptage (gènes × échantillons)
# ────────────────────────────────────────────────────────────
mat_counts <- as.matrix(round(cts_df[ , -1]))
rownames(mat_counts) <- cts_df$gene_id

# ────────────────────────────────────────────────────────────
# Objet DESeq2 et contraste pairwise
# ────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = mat_counts,
                              colData   = meta,
                              design    = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", cond1, cond2)) %>%
  as_tibble(rownames = "gene_id") %>%
  arrange(padj)

# ────────────────────────────────────────────────────────────
# Écriture du TSV résultat
# ────────────────────────────────────────────────────────────
dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
write_tsv(res, out_tsv)
