#!/usr/bin/env Rscript
# =====================================================================
# gckr_phewas.R -- PheWAS dot plot for the GCKR (2p23) locus alone (the left
#   panel of "PheWAS of the two most pleiotropic loci"), redrawn at the
#   study-wide threshold (P<5e-8/81 = 6.17e-10) only.
#
# For each trait with >=1 GCKR lead SNP at P<5e-8, takes the strongest
# (max -log10 P) association; the plot keeps only traits that also clear the
# study-wide threshold. Y-axis trait labels and dots are both colored by
# clinical category, with a category legend (no more filled/open distinction,
# since every plotted point is study-wide significant by construction).
#
# Input : data/ST2.csv   all lead SNPs at P<5e-8 (from build_study_wide_tables.R's
#                             export of workbook sheet ST2)
# Output: out/GCKR_phewas.csv          per-trait table, ALL 33 genome-wide traits
#                                      (StudyWide column flags the 25 that are plotted)
#         out/GCKR_phewas_readme.tsv
#         out/GCKR_phewas.png          the plot (study-wide traits only)
#
# Run   : Rscript scripts/gckr_phewas.R
#
# Base R only -- no CRAN packages reachable from this environment.
# =====================================================================

Sys.setlocale("LC_COLLATE", "C")

P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS

IN_ST2  <- "data/ST2.csv"
OUT_CSV <- "out/GCKR_phewas.csv"
OUT_DOC <- "out/GCKR_phewas_readme.tsv"
OUT_PNG <- "out/GCKR_phewas.png"

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

# ---- load & compute -------------------------------------------------------
st2 <- read.csv(IN_ST2, stringsAsFactors = FALSE, check.names = FALSE)
st2$neglog10P <- neglog10p(st2$BETA, st2$SE)

gckr <- st2[st2$NearestGene == "GCKR", ]
message("GCKR lead-SNP rows at P<5e-8: ", nrow(gckr), " (", length(unique(gckr$Trait)), " traits)")

# one row per trait: the strongest (max -log10 P) association at the locus
per_trait <- do.call(rbind, lapply(split(gckr, gckr$Trait), function(g) {
  g[which.max(g$neglog10P), c("Trait", "Category", "ID", "BETA", "SE", "P", "neglog10P")]
}))
per_trait$StudyWide <- per_trait$neglog10P > -log10(P_STUDYWIDE)
per_trait$Top_P <- ifelse(per_trait$P > 0, per_trait$P, NA_real_)   # blank where P underflowed

stopifnot(nrow(per_trait) == 33, sum(per_trait$Category == "Metabolism") > 0)

# Category display order: matches the published two-panel figure -- not
# alphabetical, but the fixed order the figure already uses (largest total
# signal first). Traits are then ranked by -log10P descending within category.
CATEGORY_ORDER <- c("Metabolism", "Liver", "Kidney", "Hematology", "Coagulation",
                     "Inflammatory", "Hormone", "Anthropometric", "Electrolyte", "Protein")
missing_cats <- setdiff(unique(per_trait$Category), CATEGORY_ORDER)
if (length(missing_cats)) stop("category not in CATEGORY_ORDER: ", paste(missing_cats, collapse = ", "))

per_trait$cat_rank <- match(per_trait$Category, CATEGORY_ORDER)
per_trait <- per_trait[order(per_trait$cat_rank, -per_trait$neglog10P), ]
rownames(per_trait) <- NULL

message("categories present: ", length(unique(per_trait$Category)), " of 10")
message("pass study-wide: ", sum(per_trait$StudyWide), " of ", nrow(per_trait))

write.csv(per_trait[, c("Trait", "Category", "ID", "BETA", "SE", "Top_P", "neglog10P", "StudyWide")],
          OUT_CSV, row.names = FALSE, na = "")
message("wrote ", OUT_CSV)

write_readme(OUT_DOC, "GCKR PheWAS panel - README", list(
  "Content"    = sprintf("One row per trait with >=1 GCKR (2p23) lead SNP at P<5e-8 (%d traits, %d categories); the strongest association per trait. %d of %d also clear the study-wide threshold.",
                          nrow(per_trait), length(unique(per_trait$Category)),
                          sum(per_trait$StudyWide), nrow(per_trait)),
  "Source"     = sprintf("%s, rows with NearestGene == GCKR.", IN_ST2),
  "Method"     = "Grouped by Trait; strongest (max -log10 P, computed from BETA/SE) association kept per trait. neglog10P computed in log space so P==0 (underflow) rows still rank correctly.",
  "Thresholds" = sprintf("Genome-wide P < 5e-8 (dashed line, GWS). Study-wide P < 5e-8/%d = %s (solid line, SWS) -- filled points pass both, open points pass only GWS.",
                          N_EFF_TRAITS, signif(P_STUDYWIDE, 4)),
  "Category order" = paste(CATEGORY_ORDER, collapse = " > "),
  "Script"     = "scripts/gckr_phewas.R",
  "Generated"  = format(Sys.Date())
))
message("wrote ", OUT_DOC)

# ---- plot -------------------------------------------------------------
# Category palette -- visually close to the published two-panel figure
# (dark green Metabolism ... teal Tumor marker); GCKR only touches the first
# 10 of the workbook's 12 categories, so only those are needed here.
CATEGORY_COLOR <- c(
  Metabolism      = "#1a7a3c",
  Liver           = "#a6d96a",
  Kidney          = "#d4a017",
  Hematology      = "#7b5aa6",
  Coagulation     = "#9ecae1",
  Inflammatory    = "#e6339a",
  Hormone         = "#e6772e",
  Anthropometric  = "#3b5b8c",
  Electrolyte     = "#b6e2e0",
  Protein         = "#1f6fb2"
)
LABEL_INK <- "#1a1a1a"
GRID_INK  <- "#d9d9d9"
GWS_INK   <- "#8a8a8a"
SWS_INK   <- "#1a1a1a"

# Study-wide only: drop the genome-wide-only rows entirely.
d <- per_trait[per_trait$StudyWide, ]
d <- d[order(d$cat_rank, -d$neglog10P), ]
n <- nrow(d)
y <- seq(n, 1)   # top-to-bottom in table order
cats_present <- CATEGORY_ORDER[CATEGORY_ORDER %in% unique(d$Category)]  # keep the fixed order, drop absentees

message("plotting ", n, " study-wide traits across ", length(cats_present), " categories")

XCAP <- 65                                  # cap the axis like the published panel (TG's 407 goes off-scale)
x_plot <- pmin(d$neglog10P, XCAP - 2)
off_scale <- d$neglog10P > (XCAP - 2)
col_pt <- CATEGORY_COLOR[d$Category]

png(OUT_PNG, width = 1150, height = 1150, res = 150)
par(mar = c(4.2, 6.5, 3.6, 1.5), family = "sans")

plot(NA, xlim = c(0, XCAP), ylim = c(0.5, n + 0.5), axes = FALSE,
     xlab = expression(-log[10](italic(P))), ylab = "")
abline(v = pretty(c(0, XCAP)), col = GRID_INK, lwd = 0.7)
abline(v = -log10(P_STUDYWIDE), col = SWS_INK, lty = 1, lwd = 1.1)
text(-log10(P_STUDYWIDE), n + 1.1, "SWS", col = SWS_INK, cex = 0.75, font = 2, xpd = TRUE, adj = c(-0.2, 0))

axis(1, at = pretty(c(0, XCAP)), col.axis = LABEL_INK, cex.axis = 0.85)
# Y-axis trait labels colored by category, same encoding as the dots.
axis(2, at = y, labels = FALSE, tick = FALSE)
for (i in seq_len(n)) {
  mtext(d$Trait[i], side = 2, at = y[i], las = 1, cex = 0.78, line = 0.3,
        col = col_pt[i], font = 2)
}
box(bty = "l", col = LABEL_INK)

points(x_plot, y, pch = 19, col = col_pt, cex = 1.5)
# arrow marker for points censored at the axis cap (e.g. TG)
if (any(off_scale)) {
  arrows(x_plot[off_scale] - 0.01, y[off_scale], x_plot[off_scale] + 2.2, y[off_scale],
         length = 0.06, col = col_pt[off_scale], lwd = 1.6)
}

title(main = sprintf("GCKR (2p23) -- study-wide significant\n%d traits, %d categories (P < 5x10-8/%d = %s)",
                      n, length(cats_present), N_EFF_TRAITS, signif(P_STUDYWIDE, 3)),
      cex.main = 1.0, col.main = LABEL_INK)

legend("bottomright", legend = cats_present, pch = 19, col = CATEGORY_COLOR[cats_present],
       title = "Clinical category", title.adj = 0, bty = "n", cex = 0.72, pt.cex = 1.2,
       inset = c(0.01, 0.01))

dev.off()
message("wrote ", OUT_PNG)
