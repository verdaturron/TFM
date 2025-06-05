#!/usr/bin/env Rscript
## Génère packages/ensdb/EnsDb.CanisFamiliaris.v105.sqlite via AnnotationHub

target_dir <- "packages/ensdb"
db_file    <- file.path(target_dir, "EnsDb.CanisFamiliaris.v105.sqlite")

if (!dir.exists(target_dir))
  dir.create(target_dir, recursive = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

if (!requireNamespace("AnnotationHub", quietly = TRUE))
  BiocManager::install("AnnotationHub", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(AnnotationHub)
  library(ensembldb)
})

message("⏳ Téléchargement EnsDb via AnnotationHub…")
ah  <- AnnotationHub()
qry <- query(ah, c("Canis familiaris", "EnsDb", "v105"))
if (length(qry) == 0)
  stop("EnsDb v105 pour chien introuvable sur AnnotationHub.")

edb <- qry[[1]]              

message("✅ Base reçue. Sauvegarde SQLite…")
saveDb(edb, file = db_file)  

message("📦 SQLite créé : ", db_file)
