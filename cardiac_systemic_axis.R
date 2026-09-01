#!/usr/bin/env Rscript
# =====================================================================
# cardiac_systemic_axis.R -- "Cardiac-systemic axis" figure, rebuilt on the
#   study-wide data (P < 5e-8/81 = 6.17e-10).
#
# The Discussion argues that ventricular STRUCTURE behaves genetically as a
# systemic phenotype while CONTRACTILE FUNCTION behaves as a cardiac-intrinsic
# one. This figure quantifies that claim two ways:
#
# Panel a -- correlation level (LDSC, ST7; NOT threshold-dependent).
#   For each of the 12 echocardiographic traits, the mean genetic correlation
#   with the non-cardiac traits it significantly correlates with (FDR < 0.05),
#   out of 80 non-cardiac traits tested. Diverging scale: the sign flip between
#   structural and functional traits IS the axis.
#   KEY RESULT: the flip falls between structure and function, NOT between
#   "mass/wall thickness" and "dimension/function" as the Discussion states.
#   Chamber dimensions (LVDd +0.19, 24 partners; LVDs +0.20, 22) sit firmly on
#   the structural side -- more systemically coupled than PWs (+0.37 but only
#   18) or RWT (12) -- while only EF, FS and E/A ratio invert.
#
# Panel b -- locus level (ST2 at study-wide; threshold-dependent).
#   Non-cardiac traits per study-wide cardiac locus, split by whether the locus
#   serves any functional trait. Structure-only loci average 1.4 non-cardiac
#   traits, function-serving loci 0.3 -- the same axis, independently derived.
#
# Input : data/ST7.csv    pairwise genetic correlations (workbook sheet ST7)
#         data/ST2.csv   all lead SNPs at P<5e-8 (workbook sheet ST2)
#         data/ST1.csv    trait -> clinical category (workbook sheet ST1)
# Output: out/cardiac_systemic_axis.csv           panel a table
#         out/cardiac_systemic_axis_loci.csv      panel b table
#         out/cardiac_systemic_axis_readme.tsv
#         out/cardiac_systemic_axis.png
#
# Run   : Rscript scripts/cardiac_systemic_axis.R
#
# Base R only -- no CRAN packages reachable from this environment.
# Palette: diverging pair #2166ac / #b2182b validated with the dataviz skill's
# validate_palette.js (light mode) -- all six checks PASS, worst adjacent CVD
# dE 21.1 (protan), normal-vision dE 28.7. Polarity is also encoded by bar
# direction, so identity never rests on colour alone.
# =====================================================================

Sys.setlocale("LC_COLLATE", "C")

P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS
SWS          <- -log10(P_STUDYWIDE)

IN_ST6 <- "data/ST7.csv"
IN_ST2 <- "data/ST2.csv"
IN_CAT <- "data/ST1.csv"
OUT_A   <- "out/cardiac_systemic_axis.csv"
OUT_B   <- "out/cardiac_systemic_axis_loci.csv"
OUT_DOC <- "out/cardiac_systemic_axis_readme.tsv"
OUT_PNG <- "out/cardiac_systemic_axis.png"

# Trait classes, ordered structure -> function (the axis itself)
CLASSES <- list(
  list(label = "Mass",           traits = c("LVM", "LVMI")),
  list(label = "Wall thickness", traits = c("IVSd", "IVSs", "PWd", "PWs", "RWT")),
  list(label = "Dimension",      traits = c("LVDd", "LVDs")),
  list(label = "Contraction",    traits = c("EF", "FS")),
  list(label = "Diastolic",      traits = c("EARatio"))
)
CARDIAC   <- unlist(lapply(CLASSES, `[[`, "traits"))
FUNCTIONAL <- c("EF", "FS", "EARatio")

# The workbook's ST7 header is "Trait 1" / "Trait 2" / "FDR-adjusted P"; normalise
# to syntactic names so downstream code can use $Trait1 / $Trait2 / $FDR_P.
# ST1 spells the binary ocular trait "High Myopia"; ST2/ST3/ST7 use "HighMyopia".
# Alias it so the trait -> category lookup does not silently return NA.
harmonise_traits <- function(x) {
  x <- trimws(x)
  x[x == "High Myopia"] <- "HighMyopia"
  x
}

normalise_st7 <- function(d) {
  names(d)[names(d) == "Trait 1"]          <- "Trait1"
  names(d)[names(d) == "Trait 2"]          <- "Trait2"
  names(d)[names(d) == "FDR-adjusted P"]   <- "FDR_P"
  d
}

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

# ---- load ---------------------------------------------------------------
st6 <- normalise_st7(read.csv(IN_ST6, stringsAsFactors = FALSE, check.names = FALSE))
st2 <- read.csv(IN_ST2, stringsAsFactors = FALSE, check.names = FALSE)
st1 <- read.csv(IN_CAT, stringsAsFactors = FALSE, check.names = FALSE)
catmap <- setNames(st1$Category, harmonise_traits(st1$Trait))

st2$neglog10P <- neglog10p(st2$BETA, st2$SE)
off <- st6[st6$Trait1 != st6$Trait2, ]

# ---- panel a: systemic coupling per cardiac trait ------------------------
rows <- lapply(CARDIAC, function(tr) {
  m <- off[off$Trait1 == tr | off$Trait2 == tr, ]
  partner <- ifelse(m$Trait1 == tr, m$Trait2, m$Trait1)
  keep <- catmap[partner] != "Echocardiography"
  m <- m[keep, ]
  sig <- m[m$FDR_P < 0.05, ]
  data.frame(Trait = tr, N_systemic = nrow(m), N_significant = nrow(sig),
             Pct = 100 * nrow(sig) / nrow(m),
             Mean_rg = if (nrow(sig)) mean(sig$rg) else 0,
             stringsAsFactors = FALSE)
})
A <- do.call(rbind, rows)
A$Class <- unlist(lapply(CLASSES, function(cl) rep(cl$label, length(cl$traits))))
write.csv(A, OUT_A, row.names = FALSE)
message("wrote ", OUT_A)
for (i in seq_len(nrow(A)))
  message(sprintf("  %-8s %-20s n_sig=%2d/%d  mean rg=%+0.3f",
                  A$Trait[i], A$Class[i], A$N_significant[i], A$N_systemic[i], A$Mean_rg[i]))

# ---- panel b: non-cardiac traits per study-wide cardiac locus ------------
sw <- st2[st2$neglog10P > SWS, ]
sw_card <- sw[sw$Trait %in% CARDIAC, ]
loci <- sort(unique(sw_card$NearestGene))
B <- do.call(rbind, lapply(loci, function(g) {
  d <- sw[sw$NearestGene == g, ]
  ct <- intersect(CARDIAC, unique(d$Trait))
  nc <- sort(unique(d$Trait[d$Category != "Echocardiography"]))
  data.frame(Locus = g,
             Cardiac_traits = paste(ct, collapse = ","),
             Serves_function = any(ct %in% FUNCTIONAL),
             N_noncardiac = length(nc),
             Noncardiac_traits = paste(nc, collapse = ","),
             stringsAsFactors = FALSE)
}))
B <- B[order(B$Serves_function, -B$N_noncardiac, B$Locus), ]
rownames(B) <- NULL
write.csv(B, OUT_B, row.names = FALSE)
message("wrote ", OUT_B)

m_struct <- mean(B$N_noncardiac[!B$Serves_function])
m_func   <- mean(B$N_noncardiac[B$Serves_function])
message(sprintf("structure-only loci: n=%d mean non-cardiac=%.2f | function-serving: n=%d mean=%.2f",
                sum(!B$Serves_function), m_struct, sum(B$Serves_function), m_func))

write_readme(OUT_DOC, "Cardiac-systemic axis - README", list(
  "Panel a"  = sprintf("Mean genetic correlation of each echocardiographic trait with the non-cardiac traits it significantly correlates with (FDR<0.05), of %d non-cardiac traits tested. From ST7 (LDSC) - NOT threshold-dependent.", A$N_systemic[1]),
  "Panel a key result" = "The sign flip separates STRUCTURE from FUNCTION, not 'mass/wall thickness' from 'dimension/function'. Chamber dimensions (LVDd +0.19 across 24 partners; LVDs +0.20 across 22) sit on the structural side, ahead of PWs (18 partners) and RWT (12). Only EF, FS and E/A ratio invert. The Discussion's 'hypertrophy systemic / contraction intrinsic' claim is correct; its grouping of chamber dimensions with contractile function is not.",
  "Panel b"  = sprintf("Non-cardiac traits per study-wide cardiac locus. Structure-only loci (n=%d) average %.2f non-cardiac traits; loci serving any functional trait (n=%d) average %.2f.", sum(!B$Serves_function), m_struct, sum(B$Serves_function), m_func),
  "Panel b outliers" = "MC4R (5 non-cardiac: BMI, weight, height, TG, HDL-C) is the only study-wide LV mass locus and the most systemically pleiotropic. NRG1 (4: TSH, free T4, creatinine, eGFR) is the second, and is chamber-dimension-only - again placing dimensions on the structural side.",
  "MIR4662B" = "Serves LVMI *and* LVDd/LVDs/EF/FS - the only locus spanning both sides of the axis, and the only study-wide locus for LV mass index. Carries no non-cardiac trait.",
  "Thresholds" = sprintf("Panel b uses study-wide P < 5e-8/%d = %s. Panel a is LDSC and does not depend on the GWAS threshold.", N_EFF_TRAITS, signif(P_STUDYWIDE, 4)),
  "Script"   = "scripts/cardiac_systemic_axis.R",
  "Generated" = format(Sys.Date())
))
message("wrote ", OUT_DOC)

# ---- plot ---------------------------------------------------------------
POS_COL   <- "#b2182b"   # systemic coupling (positive rg)
NEG_COL   <- "#2166ac"   # inverse coupling  (negative rg)
BAR_COL   <- "#8a8f98"   # panel b, single series -> one flat colour, no legend
INK       <- "#1a1a1a"
MUTED     <- "#666666"
GRID      <- "#e6e6e6"
RULE      <- "#c9c9c9"

png(OUT_PNG, width = 1900, height = 1150, res = 150)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 1))

## -- panel a --------------------------------------------------------------
par(mar = c(4.4, 11.4, 4.4, 3.6), family = "sans")

n <- nrow(A)
# rows top-to-bottom in CLASSES order, with a blank slot between classes
slot <- numeric(n); cls <- A$Class; pos <- 0; prev <- NA
for (i in seq_len(n)) {
  if (!is.na(prev) && cls[i] != prev) pos <- pos + 0.85
  pos <- pos + 1; slot[i] <- pos; prev <- cls[i]
}
y <- max(slot) - slot + 1
XL <- c(-0.45, 0.50)

plot(NA, xlim = XL, ylim = c(0.2, max(y) + 0.9), axes = FALSE, xlab = "", ylab = "")
abline(v = pretty(XL, 6), col = GRID, lwd = 0.7)
abline(v = 0, col = RULE, lwd = 1.2)

# bars: 4px-equivalent rounded ends are not available in base R; use flat rects
# anchored to the zero baseline, with a surface gap between adjacent bars.
h <- 0.34
rect(0, y - h, A$Mean_rg, y + h,
     col = ifelse(A$Mean_rg >= 0, POS_COL, NEG_COL), border = "white", lwd = 0.8)

# direct labels: value + how many systemic partners it rests on
for (i in seq_len(n)) {
  lab <- sprintf("%+.2f  (%d)", A$Mean_rg[i], A$N_significant[i])
  at  <- A$Mean_rg[i]
  adj <- if (A$Mean_rg[i] >= 0) -0.08 else 1.08
  text(at, y[i], lab, adj = adj, cex = 0.66, col = INK, xpd = TRUE)
}

axis(1, at = pretty(XL, 6), col = RULE, col.axis = MUTED, cex.axis = 0.72)
mtext(expression("mean " * r[g] * " with systemic traits (FDR < 0.05)"),
      side = 1, line = 2.5, cex = 0.76, col = MUTED)
for (i in seq_len(n))
  mtext(A$Trait[i], side = 2, at = y[i], las = 1, line = 0.4, cex = 0.76, font = 2, col = INK)

# class brackets in the outer margin, well clear of the trait labels
BRK <- XL[1] - 0.295
for (cl in CLASSES) {
  idx <- which(A$Trait %in% cl$traits)
  yy <- range(y[idx])
  segments(BRK, yy[1] - 0.34, BRK, yy[2] + 0.34, col = RULE, lwd = 1.4, xpd = TRUE)
  text(BRK - 0.022, mean(yy), cl$label, srt = 90, cex = 0.66, font = 2, col = MUTED, xpd = TRUE)
}

title(main = "a  Structure couples to systemic traits; function does not",
      cex.main = 0.92, col.main = INK, adj = 0, line = 2.6)
mtext("Each trait's mean genetic correlation with the non-cardiac traits it significantly correlates with; (n) = how many of 80",
      side = 3, line = 1.2, cex = 0.66, col = MUTED, adj = 0)
# ASCII only -- the PNG device mangles arrows and other non-ASCII glyphs
text(0.50, max(y) + 0.55, "systemic >>", adj = 1, cex = 0.66, font = 2, col = POS_COL, xpd = TRUE)
text(-0.45, max(y) + 0.55, "<< inverse", adj = 0, cex = 0.66, font = 2, col = NEG_COL, xpd = TRUE)

## -- panel b --------------------------------------------------------------
par(mar = c(4.4, 10.6, 4.4, 8.6))

B2 <- B[order(B$Serves_function, -B$N_noncardiac, decreasing = c(TRUE, FALSE), method = "radix"), ]
B2 <- rbind(B[!B$Serves_function, ][order(-B$N_noncardiac[!B$Serves_function]), ],
            B[B$Serves_function, ][order(-B$N_noncardiac[B$Serves_function]), ])
nb <- nrow(B2)
slotb <- numeric(nb); prev <- NA; pos <- 0
for (i in seq_len(nb)) {
  if (!is.na(prev) && B2$Serves_function[i] != prev) pos <- pos + 0.85
  pos <- pos + 1; slotb[i] <- pos; prev <- B2$Serves_function[i]
}
yb <- max(slotb) - slotb + 1
XB <- c(0, max(B2$N_noncardiac) + 1.1)

plot(NA, xlim = XB, ylim = c(0.2, max(yb) + 0.9), axes = FALSE, xlab = "", ylab = "")
abline(v = 0:max(B2$N_noncardiac), col = GRID, lwd = 0.7)
rect(0, yb - 0.34, B2$N_noncardiac, yb + 0.34, col = BAR_COL, border = "white", lwd = 0.8)
for (i in seq_len(nb)) {
  if (B2$N_noncardiac[i] > 0)
    text(B2$N_noncardiac[i], yb[i], B2$Noncardiac_traits[i], adj = -0.09, cex = 0.6, col = MUTED, xpd = TRUE)
}
axis(1, at = 0:max(B2$N_noncardiac), col = RULE, col.axis = MUTED, cex.axis = 0.72)
mtext("non-cardiac traits at the locus (study-wide)", side = 1, line = 2.5, cex = 0.76, col = MUTED)
for (i in seq_len(nb))
  mtext(B2$Locus[i], side = 2, at = yb[i], las = 1, line = 0.4, cex = 0.7, font = 4, col = INK)

BRKB <- -2.95
for (grp in c(FALSE, TRUE)) {
  idx <- which(B2$Serves_function == grp)
  yy <- range(yb[idx])
  segments(BRKB, yy[1] - 0.34, BRKB, yy[2] + 0.34, col = RULE, lwd = 1.4, xpd = TRUE)
  lab <- if (grp) sprintf("serves function (mean %.1f)", m_func) else sprintf("structure only (mean %.1f)", m_struct)
  text(BRKB - 0.18, mean(yy), lab, srt = 90, cex = 0.62, font = 2, col = MUTED, xpd = TRUE)
}

title(main = "b  The same split at locus level", cex.main = 0.95, col.main = INK, adj = 0, line = 2.6)
mtext("Study-wide cardiac loci and the non-cardiac traits they carry",
      side = 3, line = 1.2, cex = 0.66, col = MUTED, adj = 0)

dev.off()
message("wrote ", OUT_PNG)
