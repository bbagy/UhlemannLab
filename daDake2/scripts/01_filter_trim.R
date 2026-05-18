#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

opt_list <- list(
  make_option("--dada2_dir",  type="character"),
  make_option("--fastq_dir",  type="character"),
  make_option("--project",    type="character"),
  make_option("--filt_sub",   type="character"),
  make_option("--trimleft_f", type="integer",   default=0),
  make_option("--trimleft_r", type="integer",   default=0),
  make_option("--trunclen_f", type="integer"),
  make_option("--trunclen_r", type="integer"),
  make_option("--type",       type="character"),
  make_option("--date",       type="character"),
  make_option("--failed_csv", type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

dada2_dir  <- opt$dada2_dir
fastq_dir  <- opt$fastq_dir
project    <- opt$project
filt_sub   <- opt$filt_sub
date       <- opt$date
type_run   <- opt$type

# ── Load failed samples ───────────────────────────────────────────────────────
failed_samples <- character(0)
if (file.exists(opt$failed_csv)) {
  fc <- read.csv(opt$failed_csv, stringsAsFactors=FALSE)
  if (nrow(fc) > 0) failed_samples <- fc$sample
}

# ── Discover R1/R2 ───────────────────────────────────────────────────────────
if (type_run == "illumina_ITS") {
  fnFs <- sort(list.files(fastq_dir, pattern="(_L001_R1_001|_R1_001|_R1)\\.fastq\\.gz$", full.names=TRUE))
  fnRs <- sort(list.files(fastq_dir, pattern="(_L001_R2_001|_R2_001|_R2)\\.fastq\\.gz$", full.names=TRUE))
  sample.names <- sub("(_L001_R1_001|_R1_001|_R1)\\.fastq\\.gz$", "", basename(fnFs))
} else {
  fnFs <- sort(list.files(fastq_dir, pattern="(_L001_R1_001|_R1_001|_R1)\\.fastq\\.gz$", full.names=TRUE))
  fnRs <- sort(list.files(fastq_dir, pattern="(_L001_R2_001|_R2_001|_R2)\\.fastq\\.gz$", full.names=TRUE))
  sample.names <- sub("(_L001_R1_001|_R1_001|_R1)\\.fastq\\.gz$", "", basename(fnFs))
}

# exclude failed
keep <- !(sample.names %in% failed_samples)
fnFs <- fnFs[keep]; fnRs <- fnRs[keep]; sample.names <- sample.names[keep]

if (length(fnFs) == 0) stop("[filter_trim] No passing samples after excluding failed.csv entries.")
cat(sprintf("[filter_trim] %d samples to process\n", length(fnFs)))

multiplot <- function(..., plotlist=NULL, cols=1, rows=1) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  i <- 1
  while (i < numPlots) {
    numToPlot <- min(numPlots - i + 1, cols * rows)
    layout <- matrix(seq(i, i + cols * rows - 1), ncol=cols, nrow=rows, byrow=TRUE)
    if (numToPlot == 1) {
      print(plots[[i]])
    } else {
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout))))
      for (j in i:(i + numToPlot - 1)) {
        matchidx <- as.data.frame(which(layout == j, arr.ind=TRUE))
        print(plots[[j]], vp=viewport(layout.pos.row=matchidx$row, layout.pos.col=matchidx$col))
      }
    }
    i <- i + numToPlot
  }
}

# ── Raw quality plot ──────────────────────────────────────────────────────────
raw_pdf <- file.path(dada2_dir, sprintf("%s.%s.qualityProfiles.pdf", project, date))
pdf(raw_pdf)
for (sn in sample.names) {
  plist <- list()
  fn <- fnFs[sample.names == sn]
  if (length(fn) && file.exists(fn)) plist[[1]] <- plotQualityProfile(fn)
  fn <- fnRs[sample.names == sn]
  if (length(fn) && file.exists(fn)) plist[[2]] <- plotQualityProfile(fn)
  if (length(plist)) tryCatch(print(multiplot(plotlist=plist, cols=1, rows=2)), error=function(e) NULL)
}
dev.off()
cat("[filter_trim] Raw quality profiles saved\n")

# ── filterAndTrim ─────────────────────────────────────────────────────────────
filt_path <- file.path(dada2_dir, filt_sub)
dir.create(filt_path, recursive=TRUE, showWarnings=FALSE)

if (type_run == "illumina_ITS") {
  filtFs <- file.path(filt_path, paste0(sample.names, "_1_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_2_filt.fastq.gz"))
  filter_out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen   = c(opt$trunclen_f, opt$trunclen_r),
    maxEE      = c(2, 5), truncQ=2, rm.phix=TRUE,
    maxN=0, minLen=50, compress=TRUE, multithread=TRUE,
    verbose=TRUE, matchIDs=TRUE
  )
} else {
  filtFs <- file.path(filt_path, paste0(sample.names, "_R1_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_R2_filt.fastq.gz"))
  filter_out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen  = c(opt$trunclen_f, opt$trunclen_r),
    trimLeft  = c(opt$trimleft_f, opt$trimleft_r),
    maxEE     = c(2, 2), truncQ=2, rm.phix=TRUE,
    compress  = TRUE, multithread=TRUE, verbose=TRUE, matchIDs=TRUE
  )
}

# log samples that filtered to zero reads
zero_out <- rowSums(filter_out) == 0 | filter_out[,"reads.out"] == 0
if (any(zero_out)) {
  cat(sprintf("[filter_trim] WARNING: %d sample(s) have 0 reads after filtering: %s\n",
              sum(zero_out), paste(sample.names[zero_out], collapse=", ")))
  # append to failed.csv
  extra <- data.frame(
    sample   = sample.names[zero_out],
    r1       = as.character(fnFs[zero_out]),
    r2       = as.character(fnRs[zero_out]),
    r1_bytes = NA, r2_bytes = NA,
    reason   = "zero_reads_after_filter",
    stringsAsFactors = FALSE
  )
  prev <- read.csv(opt$failed_csv, stringsAsFactors=FALSE)
  write.csv(rbind(prev, extra), opt$failed_csv, row.names=FALSE, quote=FALSE)
}

# save filter stats RDS for track.csv later
saveRDS(filter_out, file.path(dada2_dir, "2_rds", "filter_stats.rds"))
cat("[filter_trim] filter_stats.rds saved\n")

# ── Filtered quality plot ─────────────────────────────────────────────────────
filt_pdf <- file.path(dada2_dir, sprintf("%s.%s.qualityProfiles.filt.pdf", project, date))
existing_filtFs <- filtFs[file.exists(filtFs) & file.info(filtFs)$size > 100]
existing_filtRs <- filtRs[file.exists(filtRs) & file.info(filtRs)$size > 100]
pdf(filt_pdf)
for (i in seq_along(existing_filtFs)) {
  plist <- list()
  tryCatch({ plist[[1]] <- plotQualityProfile(existing_filtFs[i]) }, error=function(e) NULL)
  tryCatch({ plist[[2]] <- plotQualityProfile(existing_filtRs[i]) }, error=function(e) NULL)
  if (length(plist)) tryCatch(print(multiplot(plotlist=plist, cols=1, rows=2)), error=function(e) NULL)
}
dev.off()
cat("[filter_trim] Filtered quality profiles saved\n")

# done marker for Snakemake
writeLines("done", file.path(filt_path, ".done"))
cat("[filter_trim] COMPLETE\n")
