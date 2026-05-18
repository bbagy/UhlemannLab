#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

opt_list <- list(
  make_option("--dada2_dir",   type="character"),
  make_option("--project",     type="character"),
  make_option("--date",        type="character"),
  make_option("--db",          type="character"),
  make_option("--sdb",         type="character", default=""),
  make_option("--remove_host", type="character", default="TRUE"),
  make_option("--type",        type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

dada2_dir   <- opt$dada2_dir
project     <- opt$project
date        <- opt$date
rds_dir     <- file.path(dada2_dir, "2_rds")
remove_host <- toupper(opt$remove_host) == "TRUE"

seqtab.nochim <- readRDS(file.path(rds_dir, sprintf("seqtab.nochim.%s.%s.rds", project, date)))
cat(sprintf("[taxonomy] assignTaxonomy with %s\n", opt$db))

tax <- assignTaxonomy(
  seqtab.nochim, opt$db,
  taxLevels   = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
  minBoot     = 80,
  verbose     = FALSE,
  multithread = TRUE
)

if (nzchar(opt$sdb) && file.exists(opt$sdb)) {
  cat("[taxonomy] addSpecies with species DB\n")
  tax <- addSpecies(tax, opt$sdb)
}

# ── Remove NA phylum, Chloroplast, Mitochondria (16S only) ───────────────────
if (remove_host) {
  n_before <- nrow(tax)

  is_na_phylum  <- is.na(tax[, "Phylum"])
  tax <- tax[!is_na_phylum, , drop=FALSE]
  seqtab.nochim <- seqtab.nochim[, !is_na_phylum, drop=FALSE]

  is_chloro <- !is.na(tax[, "Order"])  & tax[, "Order"]  == "Chloroplast"
  tax <- tax[!is_chloro, , drop=FALSE]
  seqtab.nochim <- seqtab.nochim[, !is_chloro, drop=FALSE]

  is_mito  <- !is.na(tax[, "Family"]) & tax[, "Family"] == "Mitochondria"
  tax <- tax[!is_mito, , drop=FALSE]
  seqtab.nochim <- seqtab.nochim[, !is_mito, drop=FALSE]

  cat(sprintf("[taxonomy] Host removal: %d → %d ASVs\n", n_before, nrow(tax)))

  # save cleaned seqtab back (needed for export step)
  saveRDS(seqtab.nochim,
          file.path(rds_dir, sprintf("seqtab.nochim.%s.%s.rds", project, date)))
}

saveRDS(tax, file.path(rds_dir, sprintf("tax.%s.%s.rds", project, date)))
cat("[taxonomy] COMPLETE\n")
