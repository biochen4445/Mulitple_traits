#!/usr/bin/env Rscript
# =====================================================================
# build_st6_all_sig.R -- build Supplementary Table S6 (per-variant
#   pleiotropy) from data/all_sig.txt.gz, restricted to STUDY-WIDE
#   SIGNIFICANT SNPs (P < 5e-8/81 = 6.17e-10).
#
# Companion to build_st5_all_sig.R (per-nearest-gene locus). Same input, same
# threshold, same hygiene; the unit of analysis is the VARIANT instead of the
# locus.
#
# History
#   This build uses every study-wide significant variant (282,674), so a variant's trait
#   count is what that variant is actually associated with. 
#   Top SNP-level hub: rs1260326 (GCKR) 25 traits here,
#   versus rs77768175 (HECTD4) 21 in the lead-SNP build.
#
# Input : ../data/all_sig.txt.gz   ALL trait-SNP associations at P<5e-8
#                               (tab-delimited, gzipped)
#   
# Output: ../out/ST6_studywide_allSNP.csv
#         ../out/ST6_studywide_allSNP_2traits.csv
#         ../out/ST6_studywide_allSNP_readme.tsv
#         
# Run   : Rscript Script/build_st6_all_sig.R
#
# Base R only -- no CRAN packages reachable from this environment.
# =====================================================================

# String sorting must not depend on the machine's locale, or every Rank shifts
# between environments. C collation is byte order, which is what the published
# workbook and Python tooling used.
Sys.setlocale("LC_COLLATE", "C")

# ---- constants -------------------------------------------------------
P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81                              # Li & Ji (2005), null-SNP Z correlation matrix
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS     # 6.172840e-10
SWS          <- -log10(P_STUDYWIDE)

IN_SIG      <- "../data/all_sig.txt.gz"
OUT_DIR <- "../out"
OUT_CSV <- file.path(OUT_DIR, "ST6_studywide_allSNP.csv")
OUT_CSV2 <- file.path(OUT_DIR, "ST6_studywide_allSNP_2traits.csv")
OUT_DOC <- file.path(OUT_DIR, "ST6_studywide_allSNP_readme.tsv")

# Analysis runs and the workbook disagree on two trait labels; the workbook
# wins (same map as build_study_wide_tables.R / build_st5_all_sig.R).
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
need <- c("ID", "CHROM", "POS", "Ref", "A1", "A1FREQ", "BETA", "SE", "P",
          "Trait", "NearestGene", "Category")
missing_cols <- setdiff(need, names(sig))
if (length(missing_cols)) stop("missing column(s) in ", IN_SIG, ": ", paste(missing_cols, collapse = ", "))
message("read ", nrow(sig), " genome-wide significant associations from ", IN_SIG)

# ---- hygiene -----------------------------------------------------------
sig$Trait       <- trimws(sig$Trait)
sig$NearestGene <- trimws(sig$NearestGene)
sig$ID          <- trimws(sig$ID)

hit <- sig$Trait %in% names(CANON_TRAIT)
if (any(hit)) {
  tb <- table(sig$Trait[hit])
  message("canonicalised trait labels: ",
          paste(sprintf("%s->%s (%d rows)", names(tb), CANON_TRAIT[names(tb)], as.integer(tb)),
                collapse = ", "))
  sig$Trait[hit] <- CANON_TRAIT[sig$Trait[hit]]
}
hitg <- sig$NearestGene %in% names(GENE_FROM_SERIAL)
if (any(hitg)) {
  message("repaired Excel date-serial gene symbols: ", sum(hitg), " rows")
  sig$NearestGene[hitg] <- GENE_FROM_SERIAL[sig$NearestGene[hitg]]
}

stopifnot(!any(is.na(sig$BETA)), !any(is.na(sig$SE)), all(sig$SE > 0),
          !any(is.na(sig$ID)), !any(sig$ID == ""))

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
        "), ", length(unique(sig$ID)), " unique variants, ", length(unique(sig$Trait)), " traits")

# ---- collapse multi-allelic records -------------------------------------
# A multi-allelic variant appears once per alternate allele, so the same
# (SNP, trait) can occur on several rows. An "association" is one SNP-trait
# pair, so the strongest allele is kept and the rest dropped.
key <- paste(sig$ID, sig$Trait, sep = "\r")
n_multi <- sum(duplicated(key))
if (n_multi) {
  sig <- sig[order(key, -sig$neglog10P), ]
  sig <- sig[!duplicated(paste(sig$ID, sig$Trait, sep = "\r")), ]
  message("collapsed ", n_multi, " extra allele record(s) at multi-allelic sites ",
          "(strongest allele kept per SNP-trait pair)")
}
stopifnot(!any(duplicated(paste(sig$ID, sig$Trait, sep = "\r"))))

# ---- per-variant invariants ---------------------------------------------
# CHROM/POS/Ref/NearestGene must be constant within a variant; A1 is not
# (64 variants remain multi-allelic across traits) and A1FREQ is not (it is
# the per-trait cohort frequency, and the cohorts differ by trait). Both are
# therefore taken from the variant's STRONGEST association, not from an
# arbitrary row.
const_ok <- function(col) all(tapply(sig[[col]], sig$ID, function(x) length(unique(x))) == 1)
stopifnot(const_ok("CHROM"), const_ok("POS"), const_ok("Ref"), const_ok("NearestGene"))
n_multi_a1 <- sum(tapply(sig$A1, sig$ID, function(x) length(unique(x))) > 1)
n_var_freq <- sum(tapply(sig$A1FREQ, sig$ID, function(x) length(unique(x))) > 1)
message(n_multi_a1, " variant(s) carry >1 alternate allele across traits; ",
        n_var_freq, " have a trait-dependent A1FREQ. ",
        "A1/A1FREQ are reported from the strongest association.")

# ---- build one row per variant ------------------------------------------
# Schema and ranking follow build_study_wide_tables.R::build_st5(), except that
# Lead_SNP becomes SNP: this build is not restricted to clumped lead variants.
per_snp <- function(g) {
  top    <- g[which.max(g$neglog10P), ]
  traits <- sort(unique(g$Trait))
  cats   <- sort(unique(g$Category))
  data.frame(
    SNP           = g$ID[1],
    CHR           = g$CHROM[1],
    POS           = g$POS[1],
    Ref           = g$Ref[1],
    A1            = top$A1,          # from the strongest association
    A1FREQ        = top$A1FREQ,      # trait-dependent; from the strongest association
    Nearest_gene  = g$NearestGene[1],
    N_traits      = length(traits),
    N_categories  = length(cats),
    Categories    = paste(cats,   collapse = "; "),
    Traits        = paste(traits, collapse = "; "),
    Top_trait     = top$Trait,
    Top_category  = top$Category,
    Top_neglog10P = round(top$neglog10P, 2),
    Top_P         = if (top$P > 0) top$P else NA_real_,   # blank, not 0, where P underflowed
    Top_BETA      = top$BETA,
    Top_SE        = top$SE,
    Pleiotropic   = if (length(traits) >= 2) "Yes" else "No",
    stringsAsFactors = FALSE
  )
}
st6 <- do.call(rbind, lapply(split(sig, sig$ID), per_snp))
# Ranking matches the published ST6: N_traits desc, Top_neglog10P desc, SNP id asc.
st6 <- st6[order(-st6$N_traits, -st6$Top_neglog10P, st6$SNP), ]
st6 <- cbind(Rank = seq_len(nrow(st6)), st6)
rownames(st6) <- NULL

stopifnot(
  nrow(st6) == length(unique(sig$ID)),
  sum(st6$N_traits) == nrow(sig),
  !any(is.na(st6$SNP)),
  !any(duplicated(st6$SNP)),
  all(st6$N_traits >= 1),
  all(st6$N_categories <= st6$N_traits),
  all(st6$Top_neglog10P > SWS)
)

message("[ST6 study-wide, all-SNP unit] ", nrow(st6), " variants, ",
        sum(st6$Pleiotropic == "Yes"), " pleiotropic, top: ",
        paste(sprintf("%s/%s (%d)", head(st6$SNP, 4), head(st6$Nearest_gene, 4),
                      head(st6$N_traits, 4)), collapse = " > "))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(st6, OUT_CSV, row.names = FALSE, na = "")
message("wrote ", OUT_CSV, " (", round(file.size(OUT_CSV) / 1e6, 1), " MB)")

ST6_2Ttraits <- st6[st6$N_traits >= 2, ]
write.csv(ST6_2Ttraits, OUT_CSV2, row.names = FALSE, na = "")

# ---- README --------------------------------------------------------------
top10 <- head(st6, 10)
write_readme(OUT_DOC, "Supplementary Table S6 (all-SNP unit, study-wide) - README", list(
  "Content"      = sprintf("One row per variant with >=1 study-wide significant association: %d variants, %d pleiotropic (>=2 traits), max %d traits. Top: %s.",
                           nrow(st6), sum(st6$Pleiotropic == "Yes"), max(st6$N_traits),
                           paste(sprintf("%s (%s, %d traits)", top10$SNP[1:5],
                                         top10$Nearest_gene[1:5], top10$N_traits[1:5]),
                                 collapse = ", ")),
  "Source"       = sprintf("%s: 696,150 records at P<5e-8, filtered to P < %s -> %d SNP-trait associations, %d unique variants, %d traits.",
                           IN_SIG, signif(P_STUDYWIDE, 4), nrow(sig),
                           length(unique(sig$ID)), length(unique(sig$Trait))),
  "UNIT -- READ THIS" = "This table covers EVERY study-wide significant variant. The first column is SNP.",
  "LD CAVEAT -- READ THIS TOO" = "Rows are NOT independent signals: neighbouring variants in an LD block carry essentially the same association, so a strong locus contributes hundreds of near-identical rows. This table answers 'what is THIS variant associated with', which is the question rs671/ALDH2-style statements in the Results text ask.",
  "Method"       = "Grouped by variant ID; ranked N_traits desc, Top_neglog10P desc, SNP id asc. Significance and Top_neglog10P computed from BETA/SE in log space, so associations below double-precision range rank correctly; the file's own study-wide flag column is cross-checked, not trusted. Multi-allelic sites collapsed to the strongest allele per SNP-trait pair. Trait labels canonicalised (Anti_nDNA -> Anti_dsDNA, Optometry_SE -> SE).",
  "Allele columns" = sprintf("CHR, POS, Ref and Nearest_gene are constant within a variant (asserted). A1 is not: %d variant(s) carry more than one alternate allele across traits. A1FREQ is not either: %d variants have a trait-dependent frequency, because the analysed cohort differs by trait. A1 and A1FREQ are therefore taken from the variant's STRONGEST association, not from an arbitrary row.",
                           n_multi_a1, n_var_freq),
  "Columns"      = "Rank; SNP; CHR; POS; Ref; A1; A1FREQ; Nearest_gene; N_traits; N_categories; Categories; Traits; Top_trait/Top_category/Top_neglog10P/Top_P/Top_BETA/Top_SE (the variant's strongest association; Top_P blank where P underflowed to 0); Pleiotropic.",
  "Threshold"    = sprintf("Study-wide P < 5e-8/%d = %s (-log10 P = %.3f), applied as an input filter.",
                           N_EFF_TRAITS, signif(P_STUDYWIDE, 4), SWS),
  "Script"       = "Script/build_st6_all_sig.R",
  "Generated"    = format(Sys.Date())
))
message("wrote ", OUT_DOC)
