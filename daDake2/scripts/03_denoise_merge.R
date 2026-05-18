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
rds_dir   <- file.path(dada2_dir, "2_rds")
out_dir   <- file.path(dada2_dir, "1_out")
dir.create(rds_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

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

good <- file.exists(filtFs) & file.info(filtFs)$size > 100 &
        file.exists(filtRs) & file.info(filtRs)$size > 100 &
        !(snames %in% failed_samples)

filtFs <- filtFs[good]; filtRs <- filtRs[good]; snames <- snames[good]
if (length(filtFs) == 0) stop("[denoise_merge] No valid filtered samples.")
cat(sprintf("[denoise_merge] %d samples\n", length(filtFs)))

# ── Load error models ─────────────────────────────────────────────────────────
errF <- readRDS(file.path(rds_dir, "errF.rds"))
errR <- readRDS(file.path(rds_dir, "errR.rds"))

# ── Dereplicate ───────────────────────────────────────────────────────────────
derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)
names(derepFs) <- snames
names(derepRs) <- snames

# ── DADA2 denoising ───────────────────────────────────────────────────────────
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# ── Merge pairs ───────────────────────────────────────────────────────────────
if (type_run == "illumina_ITS") {
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
} else {
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
}

# ── Sequence table + chimera removal ─────────────────────────────────────────
seqtab <- makeSequenceTable(mergers)
cat("[denoise_merge] Sequence length distribution:\n")
print(table(nchar(getSequences(seqtab))))

if (type_run == "illumina_ITS") {
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
} else {
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                       minFoldParentOverAbundance=1,
                                       multithread=TRUE, verbose=FALSE)
}

cat(sprintf("[denoise_merge] ASVs retained: %d / %d\n",
            ncol(seqtab.nochim), ncol(seqtab)))

# ── Save seqtab.nochim ────────────────────────────────────────────────────────
seqtab_path <- file.path(rds_dir, sprintf("seqtab.nochim.%s.%s.rds", project, date))
saveRDS(seqtab.nochim, seqtab_path)
cat(sprintf("[denoise_merge] seqtab.nochim saved: %s\n", seqtab_path))

# ── Save seqs.fna ─────────────────────────────────────────────────────────────
seqs    <- getSequences(seqtab.nochim)
headers <- paste0(">", seqs)
fasta   <- c(rbind(headers, seqs))
fna_path <- file.path(out_dir, sprintf("%s.%s.seqs.fna", project, date))
write(fasta, fna_path)
cat(sprintf("[denoise_merge] seqs.fna saved: %s\n", fna_path))

# ── Save dada stats for track.csv ─────────────────────────────────────────────
getN <- function(x) sum(getUniques(x))
dada_stats <- list(
  dadaFs_n  = sapply(dadaFs,  getN),
  dadaRs_n  = sapply(dadaRs,  getN),
  mergers_n = sapply(mergers, getN),
  nochim_n  = rowSums(seqtab.nochim)
)
saveRDS(dada_stats, file.path(rds_dir, "dada_stats.rds"))
cat("[denoise_merge] COMPLETE\n")
