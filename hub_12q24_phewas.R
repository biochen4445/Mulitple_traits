#!/usr/bin/env Rscript
# =====================================================================
# hub_12q24_phewas.R -- "PheWAS of GCKR, ABO and the 12q24 hub", redrawn for
#   just the 12q24 block: CUX2, HECTD4, RPH3A and ALDH2 (the four genes the
#   Results text treats as one LD block -- "Within this block, ... rs671
#   (ALDH2 ...)" immediately follows the HECTD4/CUX2/RPH3A sentence), and at
#   the study-wide threshold (P < 5e-8/81 = 6.17e-10) rather than the
#   genome-wide P<5e-8 the original bar counts (28/29/22/11) were drawn at.
#
# Top panel : bar chart, N_traits per gene at the study-wide threshold
#             (locus-level, from Table S4 -- the same unit Fig. 3 uses).
# Bottom panel: one row per trait that reaches study-wide significance in
#             >=1 of the 4 genes; one dot per gene, colored by -log10(P)
#             (binned; gray = doesn't reach study-wide for that gene). Traits
#             where nothing in the block passes study-wide are dropped
#             (9 of 35 genome-wide-only traits; see README).
#
# Input : data/ST2.csv   all lead SNPs at P<5e-8 (workbook sheet ST2)
#         data/ST5.csv         per-nearest-gene table, study-wide threshold
#                              (written by build_study_wide_tables.R)
# Output: out/hub_12q24_phewas.csv
#         out/hub_12q24_phewas_readme.tsv
#         out/hub_12q24_phewas.png
#
# Run   : Rscript scripts/hub_12q24_phewas.R
#
# Base R only -- no CRAN packages reachable from this environment.
# =====================================================================

Sys.setlocale("LC_COLLATE", "C")

P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS
SWS          <- -log10(P_STUDYWIDE)

GENES <- c("CUX2", "HECTD4", "RPH3A", "ALDH2")   # column order, matches the original figure

IN_ST2  <- "data/ST2.csv"
IN_ST4  <- "data/ST5.csv"
OUT_CSV <- "out/hub_12q24_phewas.csv"
OUT_DOC <- "out/hub_12q24_phewas_readme.tsv"
OUT_PNG <- "out/hub_12q24_phewas.png"

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

# ---- load -------------------------------------------------------------
st2 <- read.csv(IN_ST2, stringsAsFactors = FALSE, check.names = FALSE)
st2$neglog10P <- neglog10p(st2$BETA, st2$SE)
st4 <- read.csv(IN_ST4, stringsAsFactors = FALSE, check.names = FALSE)

sub <- st2[st2$NearestGene %in% GENES, ]
message(nrow(sub), " lead-SNP rows at P<5e-8 across ", paste(GENES, collapse = ", "))

# one value per (trait, gene): the strongest -log10P at that gene for that trait
cell <- aggregate(neglog10P ~ Trait + Category + NearestGene, sub, max)
mat <- reshape(cell, idvar = c("Trait", "Category"), timevar = "NearestGene",
                direction = "wide")
names(mat) <- sub("^neglog10P\\.", "", names(mat))
for (g in GENES) if (!g %in% names(mat)) mat[[g]] <- NA_real_
mat <- mat[, c("Trait", "Category", GENES)]

max_over_genes <- apply(mat[, GENES], 1, function(r) max(r, na.rm = TRUE))
keep <- max_over_genes > SWS
message("traits with a genome-wide hit in the block: ", nrow(mat),
        "; keep (>=1 gene study-wide): ", sum(keep), "; drop (all n.s. at study-wide): ", sum(!keep))

d <- mat[keep, ]

# Row order: category (fixed display order), then by the row's strongest
# -log10P descending within category.
CATEGORY_ORDER <- c("Metabolism", "Liver", "Kidney", "Hematology", "Coagulation",
                     "Inflammatory", "Hormone", "Anthropometric", "Electrolyte",
                     "Protein", "Vital sign", "Echocardiography")
d$cat_rank <- match(d$Category, CATEGORY_ORDER)
d$row_max  <- apply(d[, GENES], 1, function(r) max(r, na.rm = TRUE))
d <- d[order(d$cat_rank, -d$row_max), ]
rownames(d) <- NULL

write.csv(d[, c("Trait", "Category", GENES)], OUT_CSV, row.names = FALSE, na = "")
message("wrote ", OUT_CSV)

# ---- gene-level N_traits (study-wide, locus-level, from ST5) ------------
gene_n <- setNames(st4$N_traits[match(GENES, st4$Locus)], GENES)
message("study-wide N_traits per gene: ", paste(sprintf("%s=%d", GENES, gene_n), collapse = ", "))

write_readme(OUT_DOC, "12q24 hub (CUX2, HECTD4, RPH3A, ALDH2) PheWAS panel - README", list(
  "Content"      = sprintf("Bar panel: study-wide N_traits per gene (locus-level, Table S4): %s. Dot panel: %d of %d genome-wide-hit traits reach study-wide significance in >=1 of the 4 genes.",
                            paste(sprintf("%s=%d", GENES, gene_n), collapse = ", "), sum(keep), nrow(mat)),
  "Genes"        = "CUX2, HECTD4, RPH3A, ALDH2 -- the 12q24 LD block the Results text discusses as one unit (HECTD4/CUX2/RPH3A pleiotropy, then rs671/ALDH2 'within this block').",
  "Dropped traits" = paste(sprintf("%s (%s)", mat$Trait[!keep], mat$Category[!keep]), collapse = "; "),
  "Cell value"   = "Strongest (max -log10 P) lead-SNP association for that trait at that gene, from ST2 (P<5e-8). NA / gray = no lead SNP at that gene for that trait, or below the study-wide threshold.",
  "Thresholds"   = sprintf("Study-wide P < 5e-8/%d = %s (-log10P = %.2f) is the gray/colored cutoff. Genome-wide P<5e-8 only sets which trait-gene pairs exist as candidate cells at all.",
                            N_EFF_TRAITS, signif(P_STUDYWIDE, 4), SWS),
  "Note"         = "ALDH2 here is locus-level (all lead SNPs nearest-gene = ALDH2, 9 traits study-wide) -- NOT the same count as 'rs671 associated with eight traits' in the Results text, which is specific to the single rs671 SNP. Table S4 / Fig. 3 use the locus-level unit.",
  "Script"       = "scripts/hub_12q24_phewas.R",
  "Generated"    = format(Sys.Date())
))
message("wrote ", OUT_DOC)

# ---- plot ---------------------------------------------------------------
CATEGORY_COLOR <- c(
  Metabolism      = "#1a7a3c", Liver          = "#a6d96a", Kidney         = "#d4a017",
  Hematology      = "#7b5aa6", Coagulation    = "#9ecae1", Inflammatory   = "#e6339a",
  Hormone         = "#e6772e", Anthropometric = "#3b5b8c", Electrolyte    = "#b6e2e0",
  Protein         = "#1f6fb2", `Vital sign`   = "#d46fb3", Echocardiography = "#555555"
)
BAR_COLOR  <- "#b7c98e"
LABEL_INK  <- "#1a1a1a"
GRID_INK   <- "#e6e6e6"
NS_COLOR   <- "#d9d9d9"
LINE_COLOR <- "#bdbdbd"

# -log10P color bins (n.s. = below study-wide; SWS is the significance floor here,
# not the conventional genome-wide 7.3, since every plotted trait already cleared
# genome-wide significance in at least one column).
BREAKS <- c(SWS, 15, 25, 50, Inf)
BIN_COLOR <- c("#5aa0d8", "#3fa34d", "#e08a2b", "#c0392b")
BIN_LABEL <- c(sprintf("%.1f-15", SWS), "15-25", "25-50", ">=50")

bin_color <- function(v) {
  out <- rep(NS_COLOR, length(v))
  for (i in seq_along(v)) {
    if (is.na(v[i]) || v[i] <= SWS) next
    b <- findInterval(v[i], BREAKS, rightmost.closed = FALSE)
    out[i] <- BIN_COLOR[min(b + 1, length(BIN_COLOR))]
  }
  out
}

n <- nrow(d)
ng <- length(GENES)
xg <- seq_len(ng)             # gene column x-positions, shared between both panels
y  <- seq(n, 1)

png(OUT_PNG, width = 1000, height = 1450, res = 150)
layout(matrix(c(1, 3, 2, 3), nrow = 2, byrow = TRUE), heights = c(1, 5), widths = c(4, 1.3))

# -- panel 1: bar chart of study-wide N_traits per gene --------------------
par(mar = c(2.6, 6.5, 3.5, 1.5), family = "sans")
plot(NA, xlim = c(0.5, ng + 0.5), ylim = c(0, max(gene_n) * 1.25), axes = FALSE,
     xlab = "", ylab = "No. of\nassociations")
axis(2, las = 1, cex.axis = 0.75, col.axis = LABEL_INK)
rect(xg - 0.32, 0, xg + 0.32, gene_n, col = BAR_COLOR, border = NA)
text(xg, gene_n, gene_n, pos = 3, cex = 0.85, col = LABEL_INK, font = 2, xpd = TRUE)
mtext(GENES, side = 1, at = xg, line = 0.9, cex = 0.85, font = 4, col = LABEL_INK)
title(main = "PheWAS of the 12q24 hub (study-wide threshold)", cex.main = 1.05,
      col.main = LABEL_INK, adj = 0.05, line = 1.8)

# -- panel 2: dot matrix ----------------------------------------------------
par(mar = c(0.6, 6.5, 0.3, 1.5))
plot(NA, xlim = c(0.5, ng + 0.5), ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "", ylab = "")
abline(v = xg, col = GRID_INK, lwd = 0.7)

col_pt <- col_hex_matrix <- sapply(GENES, function(g) bin_color(d[[g]]))

# connecting line across the row's colored (non-n.s.) cells
for (i in seq_len(n)) {
  colored <- which(col_pt[i, ] != NS_COLOR)
  if (length(colored) >= 2) {
    lines(range(xg[colored]), c(y[i], y[i]), col = LINE_COLOR, lwd = 1)
  }
}
for (j in seq_len(ng)) points(rep(xg[j], n), y, pch = 21, bg = col_pt[, j], col = "white", cex = 1.7, lwd = 0.6)

for (i in seq_len(n)) {
  mtext(d$Trait[i], side = 2, at = y[i], las = 1, cex = 0.68, line = 0.3,
        col = CATEGORY_COLOR[d$Category[i]], font = 2)
}

# -- panel 3: legend --------------------------------------------------------
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
legend(0, 0.62, legend = c("n.s.", BIN_LABEL), pch = 21, pt.bg = c(NS_COLOR, BIN_COLOR),
       col = "white", pt.cex = 1.6, cex = 0.75, bty = "n", title = expression(-log[10](italic(P))),
       title.adj = 0, y.intersp = 1.4)
legend(0, 0.18, legend = names(CATEGORY_COLOR)[names(CATEGORY_COLOR) %in% unique(d$Category)],
       text.col = CATEGORY_COLOR[names(CATEGORY_COLOR) %in% unique(d$Category)],
       bty = "n", cex = 0.68, title = "Category (row label)", title.adj = 0, y.intersp = 1.25)

dev.off()
message("wrote ", OUT_PNG)
