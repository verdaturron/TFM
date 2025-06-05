#!/usr/bin/env Rscript
## scripts/generate_DE_report.R
## Ce script produit un rapport HTML interactif (ReportingTools)
## pour vos résultats DESeq2 (avec mini-boxplots par gène).

# 0. Installer / charger les packages nécessaires

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Installer ReportingTools et DESeq2 si besoin
for (pkg in c("ReportingTools", "DESeq2", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(ReportingTools)
library(DESeq2)
library(readr)
library(dplyr)
library(tidyr)

# ───────────────────────────────────────────────────────────────────────
# 1. Charger les résultats DESeq2 et l’objet DESeqDataSet (RDS)
# ───────────────────────────────────────────────────────────────────────
# a) Lecture du TSV de résultats
res <- read_tsv("results/deseq2/deseq2_results.tsv", show_col_types = FALSE)

# b) Chargement du DESeqDataSet préalablement sauvé
#    (votre pipeline Snakemake doit avoir créé "results/deseq2/deseq2.rds")
dds_file <- "results/deseq2/deseq2.rds"
if (!file.exists(dds_file)) {
  stop("Impossible de trouver l’objet DESeq2 à : ", dds_file)
}
dds <- readRDS(dds_file)
if (!inherits(dds, "DESeqDataSet")) {
  stop("L’objet chargé n'est pas un DESeqDataSet valide.")
}

# ───────────────────────────────────────────────────────────────────────
# 2. Préparer les données pour les tableaux + boxplots
# ───────────────────────────────────────────────────────────────────────
# 2.1 Extraire les counts normalisés
normalized_counts <- counts(dds, normalized = TRUE)

# 2.2 Construire un data.frame "long" : gene_id | sample | expr_norm | condition
sample_info <- as.data.frame(colData(dds)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  select(sample, condition)  # modifiez "condition" par le nom exact de votre colonne dans colData

df_counts_long <- as.data.frame(normalized_counts) %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "expr_norm"
  ) %>%
  left_join(sample_info, by = "sample")

# 2.3 Construire le data.frame des statuts DE (log2FC, pvalue, padj)
df_DE <- res %>%
  select(gene_id, log2FoldChange, pvalue, padj)

# ───────────────────────────────────────────────────────────────────────
# 3. Créer le rapport HTML via ReportingTools
# ───────────────────────────────────────────────────────────────────────
out_dir <- "report/DE_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Initialisation du rapport
rpt <- HTMLReport(
  shortName      = "DE_results",
  reportDirectory = out_dir,
  title          = "border vs central"
)

# Publication dans le rapport
publish(
  df_DE,
  rpt,
  name        = "Border_vs_Central",
  FeatureName = "gene_id",   # colonne ID
  plotFUN     = function(g) {
    gene <- g[[1]]
    dsub <- df_counts_long %>% filter(gene_id == gene)
    # Petits boxplots par condition
    boxplot(
      expr_norm ~ condition, 
      data    = dsub,
      outline = FALSE,
      main    = "",
      xlab    = "",
      ylab    = "",
      col     = c("skyblue", "salmon"),
      border  = "gray30"
    )
  },
  column = "Image",                          # nom de la colonne d'images
  values = c("log2FoldChange", "pvalue", "padj")  # colonnes à afficher
)

# Générer le rapport final
finish(rpt)
cat("✅ Rapport DE HTML généré dans ", out_dir, "\n")
