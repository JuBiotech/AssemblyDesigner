#!/usr/bin/env Rscript
# Whole-plasmid comparison via DECIPHER::FindSynteny (seed-and-chain), as a
# PoC alternative to the monolithic pairwise2/PairwiseAligner approach that
# struggled with full-length circular plasmid comparisons (see the main
# notebook's Step 2). Coverage and identity are computed from all syntenic
# blocks found between the two input sequences -- no manual circular padding
# needed, since blocks split naturally at the assembly's arbitrary
# linearization point instead of forcing one continuous alignment path.
#
# Usage: Rscript whole_plasmid_synteny.R <design.fasta> <assembly.fasta> <output.csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript whole_plasmid_synteny.R <design.fasta> <assembly.fasta> <output.csv>")
}
design_fa <- args[1]
assembly_fa <- args[2]
out_csv <- args[3]

# BiocManager::install() defaults to a per-user library when the R installation's
# own library (under Program Files) isn't writable -- make sure that's on the search path.
user_lib <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", "4.6")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressMessages({
  library(DECIPHER)
  library(RSQLite)
})

dbConn <- dbConnect(dbDriver("SQLite"), ":memory:")
Seqs2DB(design_fa, "FASTA", dbConn, identifier = "design", verbose = FALSE)
Seqs2DB(assembly_fa, "FASTA", dbConn, identifier = "assembly", verbose = FALSE)

syn <- FindSynteny(dbConn, minScore = 30, verbose = FALSE)
# as.dist(): "Distance is defined as one minus the hit coverage for the
# shorter of the two sequences in the pair" -- i.e. 1 - coverage.
coverage <- 1 - as.numeric(as.dist(syn))[1]

aligned <- AlignSynteny(syn, dbConn, verbose = FALSE)
pair_key <- grep("design", names(aligned), value = TRUE)[1]
blocks_aln <- aligned[[pair_key]]

total_len <- 0
total_match <- 0
for (i in seq_along(blocks_aln)) {
  a <- strsplit(as.character(blocks_aln[[i]][1]), "")[[1]]
  b <- strsplit(as.character(blocks_aln[[i]][2]), "")[[1]]
  total_match <- total_match + sum(a == b)
  total_len <- total_len + length(a)
}
pid <- if (total_len > 0) 100 * total_match / total_len else NA

result <- data.frame(
  coverage_pct = round(coverage * 100, 3),
  pid_pct = round(pid, 3),
  n_blocks = length(blocks_aln),
  aligned_columns = total_len,
  matched_columns = total_match
)
write.csv(result, out_csv, row.names = FALSE)
dbDisconnect(dbConn)
