#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

opt_list <- list(
  make_option("--dada2_dir",  type="character"),
  make_option("--project",    type="character"),
  make_option("--filt_sub",   type="character"),
  make_option("--type",       type="character"),
  make_option("--date",       type="character"),
  make_option("--failed_csv", type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

dada2_dir <- opt$dada2_dir
project   <- opt$project
date      <- opt$date
type_run  <- opt$type
filt_path <- file.path(dada2_dir, opt$filt_sub)

# ── Load failed samples ───────────────────────────────────────────────────────
failed_samples <- character(0)
if (file.exists(opt$failed_csv)) {
  fc <- read.csv(opt$failed_csv, stringsAsFactors=FALSE)
  if (nrow(fc) > 0) failed_samples <- unique(fc$sample)
}

# ── Discover filtered fastqs ──────────────────────────────────────────────────
if (type_run == "illumina_ITS") {
  filtFs <- sort(list.files(filt_path, pattern="_1_filt\\.fastq\\.gz$", full.names=TRUE))
  filtRs <- sort(list.files(filt_path, pattern="_2_filt\\.fastq\\.gz$", full.names=TRUE))
  snames <- sub("_1_filt\\.fastq\\.gz$", "", basename(filtFs))
} else {
  filtFs <- sort(list.files(filt_path, pattern="_R1_filt\\.fastq\\.gz$", full.names=TRUE))
  filtRs <- sort(list.files(filt_path, pattern="_R2_filt\\.fastq\\.gz$", full.names=TRUE))
  snames <- sub("_R1_filt\\.fastq\\.gz$", "", basename(filtFs))
}

# keep only non-empty files and non-failed samples
good <- file.exists(filtFs) & file.info(filtFs)$size > 100 &
        file.exists(filtRs) & file.info(filtRs)$size > 100 &
        !(snames %in% failed_samples)

filtFs <- filtFs[good]; filtRs <- filtRs[good]; snames <- snames[good]
if (length(filtFs) == 0) stop("[learn_errors] No valid filtered samples found.")
cat(sprintf("[learn_errors] %d samples for error learning\n", length(filtFs)))

# ── learnErrors ───────────────────────────────────────────────────────────────
set.seed(100)
errF <- suppressWarnings(learnErrors(filtFs, nreads=1000000, multithread=TRUE))
set.seed(100)
errR <- suppressWarnings(learnErrors(filtRs, nreads=1000000, multithread=TRUE))

rds_dir <- file.path(dada2_dir, "2_rds")
dir.create(rds_dir, recursive=TRUE, showWarnings=FALSE)
saveRDS(errF, file.path(rds_dir, "errF.rds"))
saveRDS(errR, file.path(rds_dir, "errR.rds"))
cat("[learn_errors] errF.rds + errR.rds saved\n")

# ── Error model plots ─────────────────────────────────────────────────────────
pdf(file.path(dada2_dir, sprintf("%s.%s.splotErrors.errF1.pdf", project, date)))
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(file.path(dada2_dir, sprintf("%s.%s.splotErrors.errF2.pdf", project, date)))
plotErrors(errR, nominalQ=TRUE)
dev.off()

cat("[learn_errors] COMPLETE\n")
