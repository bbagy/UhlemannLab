#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
  library(phyloseq)
})

opt_list <- list(
  make_option("--dada2_dir", type="character"),
  make_option("--project",   type="character"),
  make_option("--date",      type="character"),
  make_option("--type",      type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

dada2_dir <- opt$dada2_dir
project   <- opt$project
date      <- opt$date
type_run  <- opt$type

rds_dir <- file.path(dada2_dir, "2_rds")
out_dir <- file.path(dada2_dir, "1_out")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# ── Load RDS ─────────────────────────────────────────────────────────────────
seqtab.nochim <- readRDS(file.path(rds_dir, sprintf("seqtab.nochim.%s.%s.rds", project, date)))
tax_raw       <- readRDS(file.path(rds_dir, sprintf("tax.%s.%s.rds", project, date)))
filter_stats  <- readRDS(file.path(rds_dir, "filter_stats.rds"))
dada_stats    <- readRDS(file.path(rds_dir, "dada_stats.rds"))

# ── track.csv ────────────────────────────────────────────────────────────────
# filterAndTrim rownames = full input paths; seqtab.nochim rownames = sample names
# → strip path and R1 suffix to align
fs_raw <- rownames(filter_stats)
fs_sn  <- basename(fs_raw)
if (type_run == "illumina_ITS") {
  fs_sn <- sub("(_L001_R1_001|_R1_001|_R1)\\.fastq\\.gz$", "", fs_sn)
} else {
  fs_sn <- sub("(_L001)?_R1(_001)?\\.fastq\\.gz$", "", fs_sn)
}
rownames(filter_stats) <- fs_sn

st_sn  <- rownames(seqtab.nochim)
fs_sub <- filter_stats[fs_sn %in% st_sn, , drop=FALSE]
fs_sub <- fs_sub[match(st_sn, rownames(fs_sub)), , drop=FALSE]

track <- cbind(
  fs_sub,
  denoisedF = dada_stats$dadaFs_n,
  denoisedR = dada_stats$dadaRs_n,
  merged    = dada_stats$mergers_n,
  nonchim   = dada_stats$nochim_n
)
colnames(track)[1:2] <- c("input", "filtered")
write.csv(track, quote=FALSE, col.names=NA,
          file=file.path(out_dir, sprintf("%s.%s.track.csv", project, date)))
cat("[export] track.csv saved\n")

# ── Merge Genus + Species (16S only) ─────────────────────────────────────────
tax <- data.frame(tax_raw, stringsAsFactors=FALSE)
if (type_run != "illumina_ITS") {
  tax$Species <- paste(tax$Genus, tax$Species)
}
colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax <- as.matrix(tax)

# ── phyloseq object ───────────────────────────────────────────────────────────
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  tax_table(tax)
)
saveRDS(ps, file.path(rds_dir, sprintf("ps.%s.%s.rds", project, date)))
cat("[export] ps.rds saved\n")

# ── CSV outputs ───────────────────────────────────────────────────────────────
otu <- t(seqtab.nochim)

write.csv(otu, quote=FALSE, col.names=NA,
          file=file.path(out_dir, sprintf("%s.%s.asv.csv", project, date)))

write.csv(tax, quote=FALSE, col.names=NA,
          file=file.path(out_dir, sprintf("%s.%s.tax.csv", project, date)))

otuTable <- cbind(otu, tax)
write.csv(otuTable, quote=FALSE, col.names=NA,
          file=file.path(out_dir, sprintf("%s.%s.asvTable.csv", project, date)))
cat("[export] asv.csv, tax.csv, asvTable.csv saved\n")

# ── Mapping files (16S only) ──────────────────────────────────────────────────
if (type_run != "illumina_ITS") {
  map_dir <- file.path(dada2_dir, "3_map")
  dir.create(map_dir, recursive=TRUE, showWarnings=FALSE)
  SampleID <- sample_names(ps)

  make_map <- function(cols, suffix) {
    df <- data.frame(matrix(ncol=length(cols), nrow=length(SampleID)),
                     stringsAsFactors=FALSE)
    colnames(df) <- cols
    df$SampleID  <- SampleID
    df[is.na(df)] <- ""
    fpath <- file.path(map_dir, sprintf("empty.%s.%s.%s.csv", date, project, suffix))
    write.csv(df, quote=FALSE, row.names=FALSE, file=fpath)
    cat(sprintf("[export] %s saved\n", basename(fpath)))
  }

  make_map(c("SampleID","StudyID","TreatmentGroup","Timepoint","Description","etc"),
           "mapping")
  make_map(c("SampleID","is_control","sample_type","sample_well"),
           "mapping.SCRub")
}

cat("[export] COMPLETE\n")
