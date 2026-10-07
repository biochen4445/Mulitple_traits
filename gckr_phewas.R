#!/usr/bin/env Rscript
# =====================================================================
# gckr_phewas.R -- PheWAS dot plot for the GCKR (2p23) locus alone (the left
#   panel of "PheWAS of the two most pleiotropic loci"), at the study-wide
#   threshold (P < 5e-8/81 = 6.17e-10) only.
#
# For each trait with >=1 study-wide significant GCKR SNP, takes the strongest
# (max -log10 P) association. Y-axis trait labels and dots are both colored by
# clinical category, with a category legend.
#
# Input : ../data/all_sig.txt.gz   ALL trait-SNP associations at P<5e-8
#                               (tab-delimited, gzipped)
#                               -- filtered here to P < 6.17e-10.
# Output: ../out/GCKR_phewas.csv          per-trait table, study-wide traits only
#         ../out/GCKR_phewas.png          the plot
#
# Run   : Rscript scripts/gckr_phewas.R
#
# Base R only -- no CRAN packages reachable from this environment.
# =====================================================================

Sys.setlocale("LC_COLLATE", "C")

P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS
SWS          <- -log10(P_STUDYWIDE)

IN_SIG  <- "../data/all_sig.txt.gz"
OUT_CSV <- "../out/GCKR_phewas.csv"
OUT_PNG <- "../out/GCKR_phewas.png"

LOCUS_GENE <- "GCKR"

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

# ---- load & filter to study-wide significant SNPs -------------------------
# all_sig.txt.gz is tab-delimited and gzipped; read.delim() opens the .gz
# transparently.
sig <- read.delim(gzfile(IN_SIG), stringsAsFactors = FALSE, check.names = FALSE)
need <- c("ID", "CHROM", "POS", "BETA", "SE", "P", "Trait", "NearestGene", "Category")
missing_cols <- setdiff(need, names(sig))
if (length(missing_cols)) stop("missing column(s) in ", IN_SIG, ": ", paste(missing_cols, collapse = ", "))
message("read ", nrow(sig), " genome-wide significant associations from ", IN_SIG)

# Significance is decided on -log10(P) recomputed from BETA/SE in log space,
# NOT on the stored P column: tens of thousands of associations in this atlas
# are at or below double-precision range (P stored as 0 or 4.94e-324) and would
# otherwise be mis-ranked. The file's own flag column is used only as a check.
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
message("kept ", nrow(sig), " study-wide significant associations (P < ",
        signif(P_STUDYWIDE, 6), "), ", length(unique(sig$Trait)), " traits")

gckr <- sig[sig$NearestGene == LOCUS_GENE, ]
if (!nrow(gckr)) stop("no study-wide significant rows with NearestGene == ", LOCUS_GENE)

# sanity: the locus must be a single contiguous block on chr2 (2p23)
stopifnot(length(unique(gckr$CHROM)) == 1)
message(LOCUS_GENE, " study-wide SNP-trait rows: ", nrow(gckr),
        " (", length(unique(gckr$ID)), " unique SNPs, ",
        length(unique(gckr$Trait)), " traits)")
message("locus span: chr", unique(gckr$CHROM), ":",
        format(min(gckr$POS), scientific = FALSE), "-",
        format(max(gckr$POS), scientific = FALSE))

# one row per trait: the strongest (max -log10 P) association at the locus
per_trait <- do.call(rbind, lapply(split(gckr, gckr$Trait), function(g) {
  g[which.max(g$neglog10P), c("Trait", "Category", "ID", "CHROM", "POS",
                              "BETA", "SE", "P", "neglog10P")]
}))
per_trait$Top_P <- ifelse(per_trait$P > 0, per_trait$P, NA_real_)   # blank where P underflowed
per_trait$N_SNPs <- as.integer(table(gckr$Trait)[per_trait$Trait])  # study-wide SNPs behind each trait

stopifnot(nrow(per_trait) == length(unique(gckr$Trait)),
          all(per_trait$neglog10P > SWS),
          sum(per_trait$Category == "Metabolism") > 0)

# Category display order: matches the number of traits within catogories and -log10p. 
# Traits are then ranked by -log10P descending within category.
# Covers all 15 atlas categories so any category the association list brings in
# (e.g. Vital sign) has a defined slot and color.
CATEGORY_ORDER <- c("Metabolism", "Hematology", "Kidney", "Liver", "Anthropometric", "Coagulation",
                    "Inflammatory", "Hormone", "Protein", "Electrolyte",
                    "Vital sign", "Cardiac marker", "Echocardiography", "Ophthalmology",
                    "Tumor marker")
missing_cats <- setdiff(unique(per_trait$Category), CATEGORY_ORDER)
if (length(missing_cats)) stop("category not in CATEGORY_ORDER: ", paste(missing_cats, collapse = ", "))

per_trait$cat_rank <- match(per_trait$Category, CATEGORY_ORDER)
per_trait <- per_trait[order(per_trait$cat_rank, -per_trait$neglog10P), ]
rownames(per_trait) <- NULL

message("categories present: ", length(unique(per_trait$Category)), " of ", length(CATEGORY_ORDER))

write.csv(per_trait[, c("Trait", "Category", "ID", "CHROM", "POS",
                        "BETA", "SE", "Top_P", "neglog10P", "N_SNPs")],
          OUT_CSV, row.names = FALSE, na = "")
message("wrote ", OUT_CSV)


# ---- plot -------------------------------------------------------------
# Category palette -- visually close to the published two-panel figure
# (dark green Metabolism ... teal Tumor marker).
CATEGORY_COLOR <- c(
  Metabolism       = "#1a7a3c",
  Liver            = "#a6d96a",
  Kidney           = "#d4a017",
  Hematology       = "#7b5aa6",
  Coagulation      = "#9ecae1",
  Inflammatory     = "#e6339a",
  Hormone          = "#e6772e",
  Anthropometric   = "#3b5b8c",
  Electrolyte      = "#b6e2e0",
  Protein          = "#1f6fb2",
  "Vital sign"     = "#c0392b",
  "Cardiac marker" = "#8c564b",
  Echocardiography = "#6b6b6b",
  Ophthalmology    = "#b07aa1",
  "Tumor marker"   = "#17a2a2"
)
stopifnot(all(CATEGORY_ORDER %in% names(CATEGORY_COLOR)))
LABEL_INK <- "#1a1a1a"
GRID_INK  <- "#d9d9d9"
SWS_INK   <- "#1a1a1a"

d <- per_trait
n <- nrow(d)
y <- seq(n, 1)   # top-to-bottom in table order
cats_present <- CATEGORY_ORDER[CATEGORY_ORDER %in% unique(d$Category)]  # keep the fixed order, drop absentees

message("plotting ", n, " study-wide traits across ", length(cats_present), " categories")

# Axis cap: keep the published panel's 65 unless the second-strongest trait
# would itself be censored, in which case widen so that only the single
# runaway signal (TG) is drawn off-scale.
XCAP <- max(65, ceiling((sort(d$neglog10P, decreasing = TRUE)[2] + 6) / 5) * 5)
x_plot <- pmin(d$neglog10P, XCAP - 2)
off_scale <- d$neglog10P > (XCAP - 2)
col_pt <- CATEGORY_COLOR[d$Category]

png(OUT_PNG, width = 1150, height = max(1150, 40 * n + 250), res = 150)
par(mar = c(4.2, 6.5, 3.6, 1.5), family = "sans")

plot(NA, xlim = c(0, XCAP), ylim = c(0.5, n + 0.5), axes = FALSE,
     xlab = expression(-log[10](italic(P))), ylab = "")
abline(v = pretty(c(0, XCAP)), col = GRID_INK, lwd = 0.7)
# The study-wide threshold is the input filter, so it is the floor of the data,
# not a cutoff separating plotted points. Drawn as a light reference marker.
abline(v = SWS, col = SWS_INK, lty = 3, lwd = 1)
text(SWS, n + 1.1, "SWS", col = SWS_INK, cex = 0.75, font = 2, xpd = TRUE, adj = c(-0.2, 0))

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

title(main = sprintf("GCKR (2p23) -- study-wide significant SNPs\n%d traits, %d categories (P < 5x10-8/%d = %s)",
                     n, length(cats_present), N_EFF_TRAITS, signif(P_STUDYWIDE, 3)),
      cex.main = 1.0, col.main = LABEL_INK)

legend("bottomright", legend = cats_present, pch = 19, col = CATEGORY_COLOR[cats_present],
       title = "Category", title.adj = 0, bty = "n", cex = 0.72, pt.cex = 1.2,
       inset = c(0.01, 0.01))

dev.off()
message("wrote ", OUT_PNG)
