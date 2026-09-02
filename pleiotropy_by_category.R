#!/usr/bin/env Rscript
# =====================================================================
# pleiotropy_by_category.R -- which clinical category multi-trait (pleiotropic)
#   study-wide signals belong to.
#
# REVISED 2026-09-02: the unit of analysis from the VARIANT,
#   and the input from the all-SNP Table S6 (data/ST6_studywide_allSNP.csv at P < 6.17e-10).
#    
#
#   *** LD CAVEAT ***
#   A variant is "pleiotropic" when >=2 distinct traits reach the study-wide
#   threshold at that variant. Variants are NOT independent: neighbouring
#   variants in an LD block carry essentially the same signal, so one strong
#   locus contributes hundreds of near-identical rows and the absolute counts
#   here are much larger than any count of independent signals. 
#
# Input : ../data/ST6_studywide_allSNP.csv   per-variant table, study-wide,
#                                         all-SNP unit (Table S6 rebuild)
#         ../data/all_sig.txt.gz             all trait-SNP associations at P<5e-8,
#                                         filtered here to the study-wide
#                                         threshold (metric 3)
# Output: ../out/pleiotropy_by_category.csv        the three-metric table
#         ../out/pleiotropy_by_category_readme.tsv
#         ../out/pleiotropy_by_category.png        bar chart (primary metric)
#
# Run   : Rscript pleiotropy_by_category.R      (cwd = this folder)
#
# Base R only (no CRAN packages reachable from this environment or, likely,
# wherever this gets rerun).
# =====================================================================

Sys.setlocale("LC_COLLATE", "C")

P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS
SWS          <- -log10(P_STUDYWIDE)

IN_ST6   <- "../data/ST6_studywide_allSNP.csv"
IN_SIG   <- "../data/all_sig.txt.gz"
OUT_CSV  <- "../out/pleiotropy_by_category.csv"
OUT_DOC  <- "../out/pleiotropy_by_category_readme.tsv"
OUT_PNG  <- "../out/pleiotropy_by_category.png"

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

# ---- load Table S6 (per variant) ----------------------------------------
st6 <- read.csv(IN_ST6, stringsAsFactors = FALSE, check.names = FALSE)
need6 <- c("SNP", "Nearest_gene", "N_traits", "N_categories", "Categories",
           "Top_trait", "Top_category", "Pleiotropic")
missing6 <- setdiff(need6, names(st6))
if (length(missing6))
  stop("missing column(s) in ", IN_ST6, ": ", paste(missing6, collapse = ", "),
       " -- is this the OLD lead-SNP ST6 (first column Lead_SNP) rather than the all-SNP rebuild?")
message("read ", nrow(st6), " study-wide variants from ", IN_ST6)

pleio <- st6[st6$Pleiotropic == "Yes", ]
stopifnot(nrow(pleio) > 0, all(pleio$N_traits >= 2))
message("pleiotropic variants (>=2 traits): ", nrow(pleio), " of ", nrow(st6),
        " (", length(unique(pleio$Nearest_gene)), " distinct nearest genes)")

# ---- metric 1: Top_category ---------------------------------------------
# The category of the single strongest (max -log10 P) association at each
# pleiotropic variant. One vote per variant -- the cleanest "how many
# multi-trait signals does each category anchor" count.
m1 <- as.data.frame(table(pleio$Top_category), stringsAsFactors = FALSE)
names(m1) <- c("Category", "Top_category_variants")

# ---- metric 2: any category touched --------------------------------------
# Every category a pleiotropic variant's traits span, counted once per variant
# per category. A variant with traits in 3 categories contributes to all 3, so
# this sums to more than the number of variants.
cats_split <- strsplit(pleio$Categories, ";\\s*")
m2 <- as.data.frame(table(unlist(cats_split)), stringsAsFactors = FALSE)
names(m2) <- c("Category", "Any_category_variants")

# ---- metric 3: trait-associations at pleiotropic variants ----------------
# Association-level count from the raw association file under the same
# study-wide filter, restricted to variants that are pleiotropic. Weights
# categories by how many of their traits actually carry a study-wide signal at
# those variants, not just whether they are present.
sig <- read.delim(gzfile(IN_SIG), stringsAsFactors = FALSE, check.names = FALSE)
need <- c("ID", "BETA", "SE", "Trait", "Category")
missing_cols <- setdiff(need, names(sig))
if (length(missing_cols)) stop("missing column(s) in ", IN_SIG, ": ", paste(missing_cols, collapse = ", "))
sig$neglog10P <- neglog10p(sig$BETA, sig$SE)          # recomputed, not the stored P
sig <- sig[sig$neglog10P > SWS, ]
# collapse multi-allelic records, exactly as Table S6 does, so the two agree
sig <- sig[!duplicated(paste(sig$ID, sig$Trait, sep = "\r")), ]
message("study-wide associations: ", nrow(sig), " (", length(unique(sig$ID)), " variants)")
stopifnot(nrow(sig) == sum(st6$N_traits), length(unique(sig$ID)) == nrow(st6))

assoc_pleio <- sig[sig$ID %in% pleio$SNP, ]
stopifnot(nrow(assoc_pleio) == sum(pleio$N_traits))
m3 <- as.data.frame(table(assoc_pleio$Category), stringsAsFactors = FALSE)
names(m3) <- c("Category", "Associations_at_pleiotropic_variants")

# ---- merge & rank ---------------------------------------------------------
out <- merge(merge(m1, m2, all = TRUE), m3, all = TRUE)
out[is.na(out)] <- 0
out <- out[order(-out$Top_category_variants, out$Category), ]
rownames(out) <- NULL

stopifnot(
  sum(out$Top_category_variants) == nrow(pleio),
  sum(out$Associations_at_pleiotropic_variants) == nrow(assoc_pleio)
)

message("ranked by Top_category_variants:")
for (i in seq_len(nrow(out))) {
  message(sprintf("  %2d. %-16s Top_category=%6d  any-category=%6d  associations=%7d",
                  i, out$Category[i], out$Top_category_variants[i],
                  out$Any_category_variants[i], out$Associations_at_pleiotropic_variants[i]))
}

# Do the three metrics agree on the ranking? The figure's claim is about order,
# so disagreement between metrics is the thing to report, not hide.
rk <- function(v) rank(-v, ties.method = "min")
ranks <- data.frame(Category = out$Category, r1 = rk(out$Top_category_variants),
                    r2 = rk(out$Any_category_variants),
                    r3 = rk(out$Associations_at_pleiotropic_variants))
disagree <- ranks$Category[ranks$r1 != ranks$r2 | ranks$r1 != ranks$r3]
message("categories whose rank differs between metrics: ",
        if (length(disagree)) paste(disagree, collapse = ", ") else "none")

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(out, OUT_CSV, row.names = FALSE)
message("wrote ", OUT_CSV)

top4 <- head(out$Category, 4)
write_readme(OUT_DOC, "Multi-trait (pleiotropic) study-wide variants by clinical category - README", list(
  "Question"   = "Does 'multi-trait signals were most frequent among hematologic traits, followed by metabolism, liver and kidney' hold at the study-wide threshold?",
  "Answer"     = sprintf("Order by the primary metric: %s. %s",
                         paste(sprintf("%d. %s (%d)", seq_len(min(6, nrow(out))),
                                       head(out$Category, 6), head(out$Top_category_variants, 6)),
                               collapse = "; "),
                         if (identical(top4, c("Hematology", "Metabolism", "Liver", "Kidney")))
                           "This matches the stated order."
                         else
                           sprintf("This does NOT match the stated order (Hematology, Metabolism, Liver, Kidney): the top four here are %s.",
                                   paste(top4, collapse = ", "))),
  "UNIT -- READ THIS" = sprintf("Unit of analysis is the VARIANT, not the locus. %d of %d study-wide variants are pleiotropic (>=2 traits), spanning %d distinct nearest genes. The earlier version of this figure counted 1,469 pleiotropic nearest-gene LOCI from the clumped lead-SNP table; the two are different units and the bar heights are not comparable.",
                         nrow(pleio), nrow(st6), length(unique(pleio$Nearest_gene))),
  "LD CAVEAT -- READ THIS TOO" = "Variants are not independent: neighbouring variants in an LD block carry essentially the same signal, so one strong locus contributes hundreds of near-identical rows. The RANKING of categories is what this figure is for; the absolute bar heights must not be quoted as numbers of loci or of independent signals.",
  "Metric agreement" = sprintf("Categories whose rank differs between the three metrics: %s.",
                         if (length(disagree)) paste(disagree, collapse = ", ") else "none"),
  "Metric 1: Top_category_variants" = "Category of the single strongest association at each pleiotropic variant. One vote per variant. Plotted.",
  "Metric 2: Any_category_variants" = "Every category a variant's traits span, counted once per variant per category (a variant can count toward >1 category, so this sums to more than the variant count).",
  "Metric 3: Associations_at_pleiotropic_variants" = "Study-wide trait-SNP associations at pleiotropic variants, by category -- weights by how many traits per category actually signal at those variants.",
  "Source"     = sprintf("%s (%d pleiotropic of %d variants) and %s (%d study-wide associations at those variants). Study-wide P < 5e-8/%d = %s, applied on -log10(P) recomputed from BETA/SE.",
                         IN_ST6, nrow(pleio), nrow(st6), IN_SIG, nrow(assoc_pleio),
                         N_EFF_TRAITS, signif(P_STUDYWIDE, 4)),
  "Script"     = "pleiotropy_by_category.R",
  "Generated"  = format(Sys.Date())
))
message("wrote ", OUT_DOC)

# ---- plot -----------------------------------------------------------------
plot_df <- out[order(out$Top_category_variants), ]  # ascending, so barplot() reads top-to-bottom descending
n <- nrow(plot_df)

# Single-series bar chart: one flat data color (no legend needed -- see
# dataviz skill, "a single series needs no legend box, the title names it").
BAR_COLOR <- "#256abf"     # sequential-blue step 500, dataviz skill default palette
LABEL_INK <- "#1a1a1a"
GRID_INK  <- "#d9d9d9"

fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

png(OUT_PNG, width = 1500, height = max(1000, 55 * n + 260), res = 150)
par(mar = c(4, 11, 5.8, 2.5), family = "sans")

xmax <- max(plot_df$Top_category_variants) * 1.18
bp <- barplot(plot_df$Top_category_variants, horiz = TRUE, names.arg = plot_df$Category,
              las = 1, col = BAR_COLOR, border = NA, xlim = c(0, xmax),
              xlab = "Pleiotropic study-wide variants (N)", cex.names = 0.95,
              cex.axis = 0.9, col.axis = LABEL_INK, col.lab = LABEL_INK,
              xaxt = "n",
              panel.first = {
                abline(v = pretty(c(0, xmax)), col = GRID_INK, lwd = 0.7)
              })
axis(1, at = pretty(c(0, xmax)), labels = fmt(pretty(c(0, xmax))),
     cex.axis = 0.9, col.axis = LABEL_INK)
text(x = plot_df$Top_category_variants, y = bp, labels = fmt(plot_df$Top_category_variants),
     pos = 4, offset = 0.4, cex = 0.85, col = LABEL_INK, xpd = TRUE)
title(main = "Multi-trait (pleiotropic) study-wide variants by clinical category",
      cex.main = 1.05, col.main = LABEL_INK, adj = 0, line = 3.4)
mtext(sprintf("Category of the strongest signal at each of %s pleiotropic variants (>=2 traits, P < 5x10-8/81)",
              fmt(nrow(pleio))),
      side = 3, line = 2.0, cex = 0.8, col = "#555555", adj = 0)

dev.off()
message("wrote ", OUT_PNG)
