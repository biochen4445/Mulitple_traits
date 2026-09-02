#!/usr/bin/env Rscript
# =====================================================================
# build_st5_all_sig.R -- build Supplementary Table S5 (per-nearest-gene
#   pleiotropy) from data/all_sig.txt.gz, restricted to STUDY-WIDE
#   SIGNIFICANT SNPs (P < 5e-8/81 = 6.17e-10).
#
# 
#   The workbook's ST5 was built one SNP per trait per locus from All_sig.txt.gz. The Fig. 4a/4b
#   PheWAS panels were rewritten on 2026-09-02 to read every study-wide
#   significant SNP. This script produces the ST5 that matches the figures.
#
#   *** This build is the ALL-SNP unit with P < 6.17e-10.***
#
#
# Input : data/all_sig.txt.gz   ALL trait-SNP associations at P<5e-8
#                               (tab-delimited, gzipped)
#        
# Output: out/ST5_studywide_allSNP.csv
#         out/ST5_studywide_allSNP_readme.tsv
#         
# Run   : Rscript scripts/build_st5_all_sig.R
#
# Base R only -- no CRAN packages reachable from this environment.
# =====================================================================

# String sorting must not depend on the machine's locale, or e.g. "Alb" and
# "ALP" swap places between environments and every Rank shifts. C collation is
# byte order, which is what the published workbook and Python tooling used.
Sys.setlocale("LC_COLLATE", "C")

# ---- constants -------------------------------------------------------
P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81                              # Li & Ji (2005), null-SNP Z correlation matrix
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS     # 6.172840e-10
SWS          <- -log10(P_STUDYWIDE)

IN_SIG  <- "data/all_sig.txt.gz"
OUT_DIR <- "out"
OUT_CSV <- file.path(OUT_DIR, "ST5_studywide_allSNP.csv")
OUT_DOC <- file.path(OUT_DIR, "ST5_studywide_allSNP_readme.tsv")

# Analysis runs and the workbook disagree on two trait labels; the workbook
# wins (same map as build_study_wide_tables.R).
CANON_TRAIT <- c(Anti_nDNA = "Anti_dsDNA", Optometry_SE = "SE", RE = "SE")

# Excel silently turned these gene symbols into date serials when a source CSV
# was opened in the GUI. Guard kept even though this input is clean.
GENE_FROM_SERIAL <- c("46082" = "MARC1", "46083" = "MARC2", "46357" = "DEC1")

# ---- helpers -----------------------------------------------------------

# Tens of thousands of P values in this atlas underflow to exactly 0 (or the
# 4.94e-324 denormal) in double precision, so P cannot be ranked or logged
# directly. log.p = TRUE keeps the tail probability in log space.
neglog10p <- function(beta, se) {
  z <- abs(beta) / se
  -(pnorm(z, lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)
}

write_readme <- function(path, title, fields) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w"); on.exit(close(con))
  writeLines(title, con)
  for (k in names(fields)) writeLines(paste(k, fields[[k]], sep = "\t"), con)
  invisible(path)
}

# ---- load --------------------------------------------------------------
sig <- read.delim(gzfile(IN_SIG), stringsAsFactors = FALSE, check.names = FALSE)
need <- c("ID", "CHROM", "POS", "A1FREQ", "BETA", "SE", "P", "Trait", "NearestGene", "Category")
missing_cols <- setdiff(need, names(sig))
if (length(missing_cols)) stop("missing column(s) in ", IN_SIG, ": ", paste(missing_cols, collapse = ", "))
message("read ", nrow(sig), " genome-wide significant associations from ", IN_SIG)

# ---- hygiene -----------------------------------------------------------
sig$Trait       <- trimws(sig$Trait)
sig$NearestGene <- trimws(sig$NearestGene)

hit <- sig$Trait %in% names(CANON_TRAIT)
if (any(hit)) {
  message("canonicalised trait labels: ",
          paste(sprintf("%s->%s (%d rows)", names(table(sig$Trait[hit])),
                        CANON_TRAIT[names(table(sig$Trait[hit]))],
                        as.integer(table(sig$Trait[hit]))), collapse = ", "))
  sig$Trait[hit] <- CANON_TRAIT[sig$Trait[hit]]
}
hitg <- sig$NearestGene %in% names(GENE_FROM_SERIAL)
if (any(hitg)) {
  message("repaired Excel date-serial gene symbols: ", sum(hitg), " rows")
  sig$NearestGene[hitg] <- GENE_FROM_SERIAL[sig$NearestGene[hitg]]
}

stopifnot(!any(is.na(sig$BETA)), !any(is.na(sig$SE)), all(sig$SE > 0),
          !any(is.na(sig$NearestGene)), !any(sig$NearestGene == ""))

# ---- filter to study-wide significance ----------------------------------
# Significance is decided on -log10(P) recomputed from BETA/SE, not on the
# stored P column. The file's own flag column is cross-checked, not trusted.
sig$neglog10P <- neglog10p(sig$BETA, sig$SE)
is_sws <- sig$neglog10P > SWS

flag_col <- grep("^Study-wide significant", names(sig), value = TRUE)
if (length(flag_col) == 1) {
  disagree <- which(is_sws != (sig[[flag_col]] == "Yes"))
  if (length(disagree)) {
    message("NOTE: ", length(disagree), " row(s) disagree with the file's '", flag_col,
            "' column (boundary rounding); recomputed -log10(P) wins:")
    print(sig[disagree, c("ID", "Trait", "NearestGene", "P", "neglog10P", flag_col)])
  }
}
sig <- sig[is_sws, ]
message("kept ", nrow(sig), " study-wide significant records (P < ", signif(P_STUDYWIDE, 6),
        "), ", length(unique(sig$ID)), " unique SNPs, ", length(unique(sig$Trait)), " traits")

# ---- collapse multi-allelic records -------------------------------------
# A multi-allelic variant appears once per alternate allele, so the same
# (SNP, trait) can occur on several rows. An "association" is one SNP-trait
# pair, so the strongest allele is kept and the rest dropped -- otherwise
# N_associations double-counts multi-allelic sites.
key <- paste(sig$ID, sig$Trait, sep = "\r")
n_multi <- sum(duplicated(key))
if (n_multi) {
  ord <- order(key, -sig$neglog10P)
  sig <- sig[ord, ]
  sig <- sig[!duplicated(paste(sig$ID, sig$Trait, sep = "\r")), ]
  message("collapsed ", n_multi, " extra allele record(s) at multi-allelic sites ",
          "(strongest allele kept per SNP-trait pair)")
}
stopifnot(!any(duplicated(paste(sig$ID, sig$Trait, sep = "\r"))))

# ---- disambiguate gene symbols reused on >1 chromosome -------------------
# A handful of small-RNA symbols (e.g. RNU6-9) are annotated at more than one
# genomic location. Grouping on the symbol alone would give a locus spanning
# two chromosomes and a Start_bp/End_bp pair that means nothing, and would also
# make Locus non-unique for lookups. Those symbols get a "(chrN)" suffix.
chr_per_gene <- tapply(sig$CHROM, sig$NearestGene, function(x) length(unique(x)))
multi_chr <- names(chr_per_gene)[chr_per_gene > 1]
sig$Locus <- sig$NearestGene
if (length(multi_chr)) {
  i <- sig$NearestGene %in% multi_chr
  sig$Locus[i] <- sprintf("%s (chr%s)", sig$NearestGene[i], sig$CHROM[i])
  message("gene symbol(s) reused on >1 chromosome, disambiguated: ",
          paste(multi_chr, collapse = ", "))
}
stopifnot(all(tapply(sig$CHROM, sig$Locus, function(x) length(unique(x))) == 1))

# ---- build one row per locus --------------------------------------------
# this build counts every study-wide significant SNP.
per_gene <- function(g) {
  top    <- g[which.max(g$neglog10P), ]
  traits <- sort(unique(g$Trait))
  cats   <- sort(unique(g$Category))
  data.frame(
    Locus          = g$Locus[1],
    CHR            = g$CHROM[1],
    Start_bp       = min(g$POS),
    End_bp         = max(g$POS),
    N_traits       = length(traits),
    N_associations = nrow(g),
    N_SNPs         = length(unique(g$ID)),
    N_categories   = length(cats),
    Categories     = paste(cats,   collapse = "; "),
    Traits         = paste(traits, collapse = "; "),
    Top_SNP        = top$ID,
    Top_trait      = top$Trait,
    Top_category   = top$Category,
    Top_neglog10P  = round(top$neglog10P, 2),
    Top_P          = if (top$P > 0) top$P else NA_real_,   # blank, not 0, where P underflowed
    Top_BETA       = top$BETA,
    Top_SE         = top$SE,
    Top_A1FREQ     = top$A1FREQ,
    Pleiotropic    = if (length(traits) >= 2) "Yes" else "No",
    stringsAsFactors = FALSE
  )
}
st5 <- do.call(rbind, lapply(split(sig, sig$Locus), per_gene))
# Ranking matches: N_traits desc, N_associations desc, locus asc.
st5 <- st5[order(-st5$N_traits, -st5$N_associations, st5$Locus), ]
st5 <- cbind(Rank = seq_len(nrow(st5)), st5)
rownames(st5) <- NULL

stopifnot(
  sum(st5$N_associations) == nrow(sig),
  !any(is.na(st5$Locus)),
  !any(duplicated(st5$Locus)),
  all(st5$N_SNPs <= st5$N_associations),
  all(st5$N_traits >= 1),
  all(st5$Start_bp <= st5$End_bp),
  all(st5$Top_neglog10P > SWS)
)

message("[ST5 study-wide, all-SNP unit] ", nrow(st5), " loci, ",
        sum(st5$Pleiotropic == "Yes"), " pleiotropic, top: ",
        paste(sprintf("%s (%d)", head(st5$Locus, 6), head(st5$N_traits, 6)), collapse = " > "))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(st5, OUT_CSV, row.names = FALSE, na = "")
message("wrote ", OUT_CSV)


# ---- README --------------------------------------------------------------
write_readme(OUT_DOC, "Supplementary Table S5 (all-SNP unit, study-wide) - README", list(
  "Content"      = sprintf("One row per nearest-gene locus: %d loci, %d pleiotropic (>=2 traits), built from every study-wide significant SNP. Top loci: %s.",
                           nrow(st5), sum(st5$Pleiotropic == "Yes"),
                           paste(sprintf("%s (%d traits)", head(st5$Locus, 6), head(st5$N_traits, 6)), collapse = ", ")),
  "Source"       = sprintf("%s: %d records at P<5e-8, filtered to P < %s -> %d SNP-trait associations, %d unique SNPs, %d traits.",
                           IN_SIG, 696150L, signif(P_STUDYWIDE, 4), nrow(sig),
                           length(unique(sig$ID)), length(unique(sig$Trait))),
  "UNIT -- READ THIS" = "This table counts EVERY study-wide significant SNP assigned to a locus.",
  "LD CAVEAT -- READ THIS TOO" = "Rows here are nearest-gene bins, NOT independent signals. Inside a long-range LD block every gene inherits the whole block's traits, so one signal is re-counted as many loci: 13 of the top 20 rows are neighbouring genes in the single 12q24 block (HECTD4, NAA25, ACAD10, ALDH2, CUX2, BRAP, PTPN11, ADAM1A, TRAFD1, ATXN2, RPL6, LINC01405, RPH3A), and the same happens at 6p21/HLA, 11q12/FADS and 2p23/GCKR-SNX17. This build gives 6,828 loci and 130 loci with >=10 traits. This table answers 'is this gene's region study-wide associated with trait X', which is what the Fig. 4a/4b PheWAS panels ask.",
  "Method"       = "Grouped by nearest gene (and chromosome, for symbols annotated at more than one locus); ranked N_traits desc, N_associations desc, locus asc. Significance and Top_neglog10P computed from BETA/SE in log space, so associations below double-precision range rank correctly; the file's own study-wide flag column is cross-checked, not trusted. Multi-allelic sites are collapsed to the strongest allele per SNP-trait pair so N_associations counts SNP-trait pairs. Trait labels canonicalised (Anti_nDNA -> Anti_dsDNA, Optometry_SE -> SE).",
  "Columns"      = "Rank; Locus; CHR; Start_bp/End_bp (span of the locus's study-wide SNPs); N_traits; N_associations (SNP-trait pairs); N_SNPs (unique study-wide SNPs); N_categories; Categories; Traits; Top_SNP/Top_trait/Top_category/Top_neglog10P/Top_P/Top_BETA/Top_SE/Top_A1FREQ (the locus's strongest association; Top_P blank where P underflowed to 0); Pleiotropic.",
  "Threshold"    = sprintf("Study-wide P < 5e-8/%d = %s (-log10 P = %.3f), applied as an input filter.",
                           N_EFF_TRAITS, signif(P_STUDYWIDE, 4), SWS),
  "Script"       = "scripts/build_st5_all_sig_SWS.R",
  "Generated"    = format(Sys.Date())
))
message("wrote ", OUT_DOC)
