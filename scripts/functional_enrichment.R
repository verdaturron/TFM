#!/usr/bin/env Rscript
## scripts/functional_enrichment.R
## Pipeline GO / KEGG pour Canis familiaris, mapping via org.Cf.eg.db

# ───────────────────────────────────────────────────────────────────────
# 0. Dépendances (installation auto si manquantes)
# ───────────────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

needed <- c(
  "clusterProfiler", "enrichplot", "org.Cf.eg.db",
  "AnnotationDbi", "readr", "dplyr", "ggplot2"
)
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Cf.eg.db)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ───────────────────────────────────────────────────────────────────────
# 1. Chemins fournis par Snakemake
# ───────────────────────────────────────────────────────────────────────
res_file <- snakemake@input[["deseq"]]   # chemin vers results/deseq2/deseq2_results.tsv
out_dir  <- snakemake@output[["dir"]]    # dossier de sortie results/enrichment
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ───────────────────────────────────────────────────────────────────────
# 2. Lecture & filtre des gènes DE
# ───────────────────────────────────────────────────────────────────────
res <- read_tsv(res_file, show_col_types = FALSE)

sig <- res |>
  filter(padj < 0.05, abs(log2FoldChange) >= 1) |>
  filter(!is.na(gene_id), !is.na(log2FoldChange))

stopifnot("gene_id" %in% colnames(sig))
gene_list <- sig$gene_id
if (length(gene_list) < 4)
  stop("Pas assez de gènes DE pour l’enrichissement.")

# ───────────────────────────────────────────────────────────────────────
# 3. Mapping ENSEMBL → SYMBOL / ENTREZ via org.Cf.eg.db
# ───────────────────────────────────────────────────────────────────────
cat("⚙️  Mapping Ensembl → SYMBOL/ENTREZ via org.Cf.eg.db…\n")

map <- AnnotationDbi::select(
  org.Cf.eg.db,
  keys       = gene_list,
  keytype    = "ENSEMBL",
  columns    = c("SYMBOL", "ENTREZID")
) %>%
  rename(
    ensembl_gene_id    = ENSEMBL,
    external_gene_name = SYMBOL,
    entrezgene_id      = ENTREZID
  ) %>%
  filter(ensembl_gene_id %in% gene_list)

if (nrow(map) == 0)
  stop("❌ org.Cf.eg.db n’a renvoyé aucun mapping pour les IDs fournis.")

sig <- sig |>
  left_join(map, by = c("gene_id" = "ensembl_gene_id")) |>
  filter(!is.na(external_gene_name) | !is.na(entrezgene_id))

symbol_list <- unique(na.omit(sig$external_gene_name))
entrez_list <- unique(na.omit(sig$entrezgene_id))

# ───────────────────────────────────────────────────────────────────────
# 4. Enrichissement GO (keyType = SYMBOL)
# ───────────────────────────────────────────────────────────────────────
ego <- enrichGO(
  symbol_list,
  OrgDb         = org.Cf.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

if (!is.null(ego) && nrow(ego@result) > 0) {
  tryCatch(
    ggsave(
      file.path(out_dir, "go_BP_dotplot.png"),
      dotplot(ego, showCategory = 20) + ggtitle("GO – Biological Process"),
      width = 8, height = 6, dpi = 300
    ),
    error = function(e) message("⚠ dotplot KO : ", e$message)
  )
  tryCatch(
    ggsave(
      file.path(out_dir, "go_BP_cnetplot.png"),
      cnetplot(ego, showCategory = 10) + ggtitle("GO BP – Cnetplot"),
      width = 8, height = 6, dpi = 300
    ),
    error = function(e) message("⚠ cnetplot KO : ", e$message)
  )
  write_tsv(as_tibble(ego), file.path(out_dir, "go_BP_enrichment.tsv"))
} else {
  message("⚠ Aucun résultat GO BP significatif.")
}

# ───────────────────────────────────────────────────────────────────────
# 5. Enrichissement KEGG (IDs ENTREZ)
# ───────────────────────────────────────────────────────────────────────
log_path <- file.path(out_dir, "log.txt")
log_conn <- file(log_path, open = "wt")
logmsg <- function(...) { cat(paste0(..., "\n"), file = log_conn); message(...) }

logmsg("Début KEGG | ENTREZ disponibles : ", length(entrez_list))

kegg <- NULL
if (length(entrez_list) > 0) {
  kegg <- tryCatch(
    enrichKEGG(
      gene          = entrez_list,
      organism      = "cfa",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    ),
    error = function(e) { logmsg("enrichKEGG KO : ", e$message); NULL }
  )
  if (!is.null(kegg) && nrow(kegg) > 0) {
    tryCatch(
      ggsave(
        file.path(out_dir, "kegg_dotplot.png"),
        dotplot(kegg, showCategory = 20) + ggtitle("KEGG pathways"),
        width = 8, height = 6, dpi = 300
      ),
      error = function(e) message("⚠ kegg dotplot KO : ", e$message)
    )
    write_tsv(as_tibble(kegg), file.path(out_dir, "kegg_enrichment.tsv"))
    logmsg("✅ KEGG OK (", nrow(kegg), " pathways).")
  } else {
    logmsg("⚠ Pas de KEGG significatif.")
  }
} else {
  logmsg("⚠ Pas d’ENTREZID après mapping.")
}
close(log_conn)

# ───────────────────────────────────────────────────────────────────────
# 6. Sauvegarde des objets R
# ───────────────────────────────────────────────────────────────────────
saveRDS(
  list(ego = ego, kegg = kegg),
  file = file.path(out_dir, "enrichment_objects.rds")
)

cat("✅ Enrichissement terminé.\n")
