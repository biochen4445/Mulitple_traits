#!/usr/bin/env Rscript
# =====================================================================
# hub_12q24_phewas.R -- "PheWAS of GCKR, ABO and the 12q24 hub" at the
#   study-wide threshold (P < 5e-8/81 = 6.17e-10).
#
# Top panel : bar chart, N_traits per gene, read from Table S5 (all-SNP unit).
# Bottom panel: one row per trait study-wide significant in >=1 of the 4
#             genes; one dot per gene, colored by -log10(P) (binned;
#             gray = no study-wide SNP at that gene for that trait).
#
# Input : data/all_sig.txt.gz   ALL trait-SNP associations at P<5e-8
#                               (tab-delimited, gzipped) -- filtered here to
#                               P < 6.17e-10. Drives the dot panel.
#         data/ST5_studywide_allSNP.csv   ST5 from all-SNP unit, study-wide. 
# Output: out/hub_12q24_phewas.csv
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

GENES <- c("CUX2", "ACAD10", "ALDH2", "NAA25", "HECTD4")   # column order, matches the original figure

IN_SIG      <- "../data/all_sig.txt.gz"                 # dot panel
IN_ST5      <- "../data/ST5_studywide_allSNP.csv"       # bar panel (all-SNP unit)

OUT_CSV <- "../out/hub_12q24_phewas.csv"

OUT_PNG <- "../out/hub_12q24_phewas.png"

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
sig <- read.delim(gzfile(IN_SIG), stringsAsFactors = FALSE, check.names = FALSE)
need <- c("ID", "CHROM", "POS", "BETA", "SE", "P", "Trait", "NearestGene", "Category")
missing_cols <- setdiff(need, names(sig))
if (length(missing_cols)) stop("missing column(s) in ", IN_SIG, ": ", paste(missing_cols, collapse = ", "))
message("read ", nrow(sig), " genome-wide significant associations from ", IN_SIG)

# Significance decided on -log10(P) recomputed from BETA/SE in log space, not on
# the stored P column (many associations are below double-precision range).
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

sub <- sig[sig$NearestGene %in% GENES, ]
if (!nrow(sub)) stop("no study-wide significant rows for ", paste(GENES, collapse = ", "))
stopifnot(length(unique(sub$CHROM)) == 1)
message(nrow(sub), " study-wide SNP-trait rows across ", paste(GENES, collapse = ", "),
        " (", length(unique(sub$ID)), " unique SNPs, ", length(unique(sub$Trait)), " traits)")
message("block span: chr", unique(sub$CHROM), ":",
        format(min(sub$POS), scientific = FALSE), "-",
        format(max(sub$POS), scientific = FALSE))

# one value per (trait, gene): the strongest -log10P at that gene for that trait
cell <- aggregate(neglog10P ~ Trait + Category + NearestGene, sub, max)
mat <- reshape(cell, idvar = c("Trait", "Category"), timevar = "NearestGene",
               direction = "wide")
names(mat) <- sub("^neglog10P\\.", "", names(mat))
for (g in GENES) if (!g %in% names(mat)) mat[[g]] <- NA_real_
mat <- mat[, c("Trait", "Category", GENES)]

# Every row is study-wide by construction now (the filter is applied to the
# input), so nothing is dropped here -- kept as an assertion instead.
row_max <- apply(mat[, GENES], 1, function(r) max(r, na.rm = TRUE))
stopifnot(all(row_max > SWS))
d <- mat
message("traits study-wide in >=1 gene of the block: ", nrow(d))

# Row order: category (fixed display order), then by the row's strongest
# -log10P descending within category.
CATEGORY_ORDER <- c("Metabolism", "Liver", "Kidney", "Hematology", "Coagulation",
                    "Inflammatory", "Hormone", "Anthropometric", "Electrolyte",
                    "Protein", "Vital sign", "Cardiac marker", "Echocardiography",
                    "Ophthalmology", "Tumor marker")
missing_cats <- setdiff(unique(d$Category), CATEGORY_ORDER)
if (length(missing_cats)) stop("category not in CATEGORY_ORDER: ", paste(missing_cats, collapse = ", "))
d$cat_rank <- match(d$Category, CATEGORY_ORDER)
d$row_max  <- apply(d[, GENES], 1, function(r) max(r, na.rm = TRUE))
d <- d[order(d$cat_rank, -d$row_max), ]
rownames(d) <- NULL

write.csv(d[, c("Trait", "Category", GENES)], OUT_CSV, row.names = FALSE, na = "")
message("wrote ", OUT_CSV)

# ---- gene-level N_traits, from Table S5 ----------------------------------
# The bar panel is a consumer of the rebuilt all-SNP Table S5, not an
# independent recount, so the table and the figure cannot drift apart.
st5 <- read.csv(IN_ST5, stringsAsFactors = FALSE, check.names = FALSE)
need_st5 <- c("Rank", "Locus", "CHR", "N_traits", "N_associations", "N_SNPs", "N_categories")
missing_st5 <- setdiff(need_st5, names(st5))
if (length(missing_st5))
  stop("missing column(s) in ", IN_ST5, ": ", paste(missing_st5, collapse = ", "),
       " -- is this the OLD lead-SNP ST5 (N_lead_SNPs/Top_lead_SNP) rather than the all-SNP build?")

i5 <- match(GENES, st5$Locus)
if (any(is.na(i5))) stop("gene(s) absent from ", IN_ST5, ": ", paste(GENES[is.na(i5)], collapse = ", "))
gene_row <- st5[i5, ]
gene_n   <- setNames(gene_row$N_traits, GENES)
message("Table S5 (all-SNP unit) per gene: ",
        paste(sprintf("%s: rank %d, %d traits / %d assoc / %d SNPs / %d categories",
                      GENES, gene_row$Rank, gene_row$N_traits, gene_row$N_associations,
                      gene_row$N_SNPs, gene_row$N_categories), collapse = "; "))

# Consistency assertion: the table must reproduce what this script's own
# filtered association list says. A mismatch means the two were built from
# different inputs or thresholds -- stop rather than plot a stale bar.
gene_n_here <- setNames(sapply(GENES, function(g) length(unique(sub$Trait[sub$NearestGene == g]))), GENES)
if (!identical(as.integer(gene_n), as.integer(gene_n_here))) {
  print(data.frame(Gene = GENES, From_TableS5 = as.integer(gene_n),
                   From_all_sig = as.integer(gene_n_here)))
  stop("Table S5 and ", IN_SIG, " disagree on N_traits -- rebuild Table S5 before redrawing.")
}
message("consistency check passed: Table S5 N_traits == recount from ", IN_SIG)


# ---- plot ---------------------------------------------------------------
CATEGORY_COLOR <- c(
  Metabolism      = "#1a7a3c", Liver          = "#a6d96a", Kidney         = "#d4a017",
  Hematology      = "#7b5aa6", Coagulation    = "#9ecae1", Inflammatory   = "#e6339a",
  Hormone         = "#e6772e", Anthropometric = "#3b5b8c", Electrolyte    = "#b6e2e0",
  Protein         = "#1f6fb2", `Vital sign`   = "#c0392b", `Cardiac marker` = "#8c564b",
  Echocardiography = "#555555", Ophthalmology = "#b07aa1", `Tumor marker` = "#17a2a2"
)
stopifnot(all(CATEGORY_ORDER %in% names(CATEGORY_COLOR)))
BAR_COLOR  <- "#b7c98e"
LABEL_INK  <- "#1a1a1a"
GRID_INK   <- "#e6e6e6"
NS_COLOR   <- "#d9d9d9"
LINE_COLOR <- "#bdbdbd"

# -log10P color bins. The study-wide threshold is the input floor, so a gray
# cell means "no study-wide SNP at this gene for this trait", not "tested and
# non-significant".
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

png(OUT_PNG, width = 1000, height = max(1450, 42 * n + 420), res = 150)
layout(matrix(c(1, 3, 2, 3), nrow = 2, byrow = TRUE), heights = c(1, 5), widths = c(4, 1.3))

# -- panel 1: bar chart of study-wide N_traits per gene --------------------
par(mar = c(2.6, 6.5, 3.5, 1.5), family = "sans")
plot(NA, xlim = c(0.5, ng + 0.5), ylim = c(0, max(gene_n) * 1.25), axes = FALSE,
     xlab = "", ylab = "No. of\ntraits")   # N_traits, not N_associations -- the old label was wrong
axis(2, las = 1, cex.axis = 0.75, col.axis = LABEL_INK)
rect(xg - 0.32, 0, xg + 0.32, gene_n, col = BAR_COLOR, border = NA)
text(xg, gene_n, gene_n, pos = 3, cex = 0.85, col = LABEL_INK, font = 2, xpd = TRUE)
mtext(GENES, side = 1, at = xg, line = 0.9, cex = 0.85, font = 4, col = LABEL_INK)
title(main = "PheWAS of the 12q24 hub (study-wide significant SNPs)", cex.main = 1.05,
      col.main = LABEL_INK, adj = 0.05, line = 1.8)

# -- panel 2: dot matrix ----------------------------------------------------
par(mar = c(0.6, 6.5, 0.3, 1.5))
plot(NA, xlim = c(0.5, ng + 0.5), ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "", ylab = "")
abline(v = xg, col = GRID_INK, lwd = 0.7)

col_pt <- sapply(GENES, function(g) bin_color(d[[g]]))

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
legend(0, 0.62, legend = c("none", BIN_LABEL), pch = 21, pt.bg = c(NS_COLOR, BIN_COLOR),
       col = "white", pt.cex = 1.6, cex = 0.75, bty = "n", title = expression(-log[10](italic(P))),
       title.adj = 0, y.intersp = 1.4)
legend(0, 0.18, legend = names(CATEGORY_COLOR)[names(CATEGORY_COLOR) %in% unique(d$Category)],
       text.col = CATEGORY_COLOR[names(CATEGORY_COLOR) %in% unique(d$Category)],
       bty = "n", cex = 0.68, title = "Category (row label)", title.adj = 0, y.intersp = 1.25)

dev.off()
message("wrote ", OUT_PNG)
