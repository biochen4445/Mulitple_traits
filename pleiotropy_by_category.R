#!/usr/bin/env Rscript
# =====================================================================
# echo_partition_rg.R -- "Echocardiographic loci partition by trait class"
#   (panel a) + "Genetic correlation with systemic markers" (panel b),
#   at the study-wide threshold only (P < 5e-8/81 = 6.17e-10).
#   Panel a is threshold-dependent (locus discovery); panel b is not --
#   genetic correlations (LDSC, ST7) do not depend on the GWAS significance
#   threshold, so panel b is unchanged by this revision.
#
# REVISION 2026-09-02: panel a's input from the full association list data/all_sig.txt.gz,
#   FILTERED TO STUDY-WIDE SIGNIFICANT SNPs before anything downstream. Every
#   significant SNP now contributes, so a locus appears whenever any study-wide SNP nearest that gene hits an echo
#   trait. Panel a grows to 33. Panel b is untouched.
#
#   *** LD CAVEAT ***
#   Rows are nearest-gene bins, not independent signals. Neighbouring genes in
#   one LD block each inherit the block's traits, so a single signal can occupy
#   several rows -- TTN and TTN-AS1 here are the same 2q31 signal under two
#   symbols. Read panel a as "which genes' regions carry study-wide echo
#   signal". 
#
# Panel a: rows = nearest-gene loci with >=1 study-wide significant SNP among
#   the 11 echocardiographic traits (LVM, LVMI, IVSd, IVSs, PWd, PWs, RWT,
#   LVDd, LVDs, EF, FS); dark cell = study-wide significant for that trait.
#   Columns continue into 6 systemic-trait categories (Adiposity, Blood
#   pressure, Metabolic, Liver, Kidney, Blood -- mapped 1:1 from the atlas's
#   Anthropometric / Vital sign / Metabolism / Liver / Kidney / Hematology);
#   shade = number of systemic traits in that category also study-wide
#   significant at the same locus (light = 1, dark = >=2).
#   Rows are grouped by primary phenotype class (more cardiac-trait hits in
#   the mass/wall-thickness block vs the dimension/function block; ties go to
#   mass). At study-wide, MC4R remains the ONLY mass/wall-thickness-primary
#   locus. Loci with hits in both blocks stay in the dimension/function block
#   but keep their mass-block cells colored, visibly bridging the two -- the
#   MIR4662B finding already flagged as a Word comment on the main text, now
#   joined by LINC00964.
#
# Panel b: rows = the 7 mass/wall-thickness traits (LVM, LVMI, IVSd, IVSs,
#   PWd, PWs, RWT) + 2 function traits (EF, FS) -- LVDd/LVDs excluded, same
#   trait set as the reference image. Columns = 10 systemic markers (ALT,
#   GGT, AST, hsCRP, ESR, uACR, uPCR, SCr, SBP, BMI). Cell = genetic
#   correlation (rg) from ST7, diverging blue-white-red; "*" = FDR < 0.05.
#
# Input : ../data/all_sig.txt.gz  ALL trait-SNP associations at P<5e-8
#                              (tab-delimited, gzipped) -- filtered here to
#                              P < 6.17e-10. Drives panel a.
#         ../data/ST7.csv         all pairwise genetic correlations (LDSC
#                              rg/SE/Z/P/FDR-adjusted P). Drives panel b;
#                              threshold-independent.
# Output: ../out/echo_partition.csv          panel a matrix (loci x traits)
#         ../out/echo_partition_readme.tsv
#         ../out/echo_rg_systemic.csv        panel b matrix (traits x markers)
#         ../out/echo_partition_rg.png       the two-panel figure
#
# Run   : Rscript scripts/echo_partition_rg.R
#
# Base R only -- no CRAN packages reachable from this environment.
# =====================================================================

Sys.setlocale("LC_COLLATE", "C")

P_GENOMEWIDE <- 5e-8
N_EFF_TRAITS <- 81
P_STUDYWIDE  <- P_GENOMEWIDE / N_EFF_TRAITS
SWS          <- -log10(P_STUDYWIDE)

CARDIAC <- c("LVM", "LVMI", "IVSd", "IVSs", "PWd", "PWs", "RWT", "LVDd", "LVDs", "EF", "FS")
MASS_TRAITS  <- c("LVM", "LVMI", "IVSd", "IVSs", "PWd", "PWs", "RWT")
DIMFN_TRAITS <- c("LVDd", "LVDs", "EF", "FS")

CAT_MAP <- c(Anthropometric = "Adiposity", `Vital sign` = "Blood pressure",
             Metabolism = "Metabolic", Liver = "Liver", Kidney = "Kidney",
             Hematology = "Blood")
SYS_COLS <- unname(CAT_MAP)   # display order

PANELB_ROWS <- c("PWs", "PWd", "IVSs", "IVSd", "LVM", "LVMI", "RWT", "EF", "FS")
PANELB_MASS <- c("PWs", "PWd", "IVSs", "IVSd", "LVM", "LVMI", "RWT")
PANELB_FUNC <- c("EF", "FS")
SYS_MARKERS <- c("ALT", "GGT", "AST", "hsCRP", "ESR", "uACR", "uPCR", "SCr", "SBP", "BMI")

IN_SIG  <- "../data/all_sig.txt.gz"
IN_ST7  <- "../data/ST7.csv"
OUT_A_CSV <- "../out/echo_partition.csv"
OUT_DOC   <- "../out/echo_partition_readme.tsv"
OUT_B_CSV <- "../out/echo_rg_systemic.csv"
OUT_PNG   <- "../out/echo_partition_rg.png"

# The workbook's ST7 header is "Trait 1" / "Trait 2" / "FDR-adjusted P"; normalise
# to syntactic names so downstream code can use $Trait1 / $Trait2 / $FDR_P.
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

# ==== Panel a: locus x trait study-wide partition ==========================
sig <- read.delim(gzfile(IN_SIG), stringsAsFactors = FALSE, check.names = FALSE)
need <- c("ID", "BETA", "SE", "P", "Trait", "NearestGene", "Category")
missing_cols <- setdiff(need, names(sig))
if (length(missing_cols)) stop("missing column(s) in ", IN_SIG, ": ", paste(missing_cols, collapse = ", "))
message("read ", nrow(sig), " genome-wide significant associations from ", IN_SIG)

# Significance is decided on -log10(P) recomputed from BETA/SE in log space,
# not the stored P column (many associations are below double-precision range).
# The file's own flag column is cross-checked, not trusted.
sig$neglog10P <- neglog10p(sig$BETA, sig$SE)
# Which echo traits have no genome-wide association AT ALL -- captured before
# the study-wide filter, so the README can tell a structurally empty column
# apart from one emptied by raising the threshold.
absent <- setdiff(CARDIAC, unique(sig$Trait))
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
# collapse multi-allelic records, as Table S5/S6 do, so counts agree with them
sig <- sig[!duplicated(paste(sig$ID, sig$Trait, sep = "\r")), ]
message("kept ", nrow(sig), " study-wide associations (", length(unique(sig$ID)), " variants, ",
        length(unique(sig$Trait)), " traits)")

sw <- sig[sig$Trait %in% CARDIAC, ]
if (!nrow(sw)) stop("no study-wide significant associations for any of the 11 echo traits")
loci <- sort(unique(sw$NearestGene))
message(length(loci), " study-wide echo loci across the 11-trait set (",
        length(unique(sw$ID)), " variants, ", nrow(sw), " associations)")
message("echo traits with NO study-wide locus: ",
        paste(setdiff(CARDIAC, unique(sw$Trait)), collapse = ", "),
        if (length(absent)) paste0("  [of which ", paste(absent, collapse = ", "),
                                   " have no genome-wide association at all]") else "")

cardiac_mat <- matrix(FALSE, nrow = length(loci), ncol = length(CARDIAC),
                      dimnames = list(loci, CARDIAC))
for (g in loci) {
  hits <- unique(sw$Trait[sw$NearestGene == g])
  cardiac_mat[g, hits] <- TRUE
}

sys_mat <- matrix(0L, nrow = length(loci), ncol = length(SYS_COLS),
                  dimnames = list(loci, SYS_COLS))
for (g in loci) {
  d <- sig[sig$NearestGene == g & sig$Category %in% names(CAT_MAP), ]
  if (nrow(d)) {
    tab <- tapply(d$Trait, CAT_MAP[d$Category], function(x) length(unique(x)))
    sys_mat[g, names(tab)] <- tab
  }
}

# primary group: more mass-block hits (ties -> mass), else dimension/function
mass_n  <- rowSums(cardiac_mat[, MASS_TRAITS,  drop = FALSE])
dimfn_n <- rowSums(cardiac_mat[, DIMFN_TRAITS, drop = FALSE])
grp <- ifelse(mass_n >= dimfn_n & mass_n > 0, "Mass / wall-thickness loci",
       ifelse(dimfn_n > 0, "Dimension / function loci", "Mass / wall-thickness loci"))

# index by name, not position -- mass_n/dimfn_n are computed before the reorder
bridging  <- names(mass_n)[mass_n > 0 & dimfn_n > 0]
mass_only <- names(mass_n)[grp[names(mass_n)] == "Mass / wall-thickness loci"]
cardiac_only <- loci[rowSums(sys_mat) == 0]

ord <- order(factor(grp, levels = c("Mass / wall-thickness loci", "Dimension / function loci")),
             -(mass_n + dimfn_n), loci)
loci   <- loci[ord]
grp    <- grp[ord]
cardiac_mat <- cardiac_mat[loci, , drop = FALSE]
sys_mat     <- sys_mat[loci, , drop = FALSE]

out_a <- data.frame(Locus = loci, Group = grp, cardiac_mat, sys_mat, check.names = FALSE)
dir.create(dirname(OUT_A_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(out_a, OUT_A_CSV, row.names = FALSE)
message("wrote ", OUT_A_CSV)

message("mass/wall-thickness-primary loci: ", paste(mass_only, collapse = ", "))
message("loci with hits in BOTH blocks (bridging): ", paste(bridging, collapse = ", "))
message("loci with NO study-wide systemic association: ", length(cardiac_only), " of ", length(loci))

describe_bridge <- function(g) sprintf("%s (mass-block: %s; dimension/function: %s)", g,
  paste(CARDIAC[cardiac_mat[g, ] & CARDIAC %in% MASS_TRAITS],  collapse = ","),
  paste(CARDIAC[cardiac_mat[g, ] & CARDIAC %in% DIMFN_TRAITS], collapse = ","))

write_readme(OUT_DOC, "Echocardiographic loci partition (study-wide, all-SNP unit) - README", list(
  "Content"     = sprintf("%d nearest-gene loci reach study-wide significance (P<5e-8/81) for >=1 of the 11 echo traits. %d are mass/wall-thickness-primary, %d are dimension/function-primary.",
                          length(loci), sum(grp == "Mass / wall-thickness loci"),
                          sum(grp == "Dimension / function loci")),
  "UNIT -- READ THIS" = sprintf("Built from EVERY study-wide significant SNP (%s), not the clumped lead SNPs. The lead-SNP version of this panel had 15 loci; this one has %d. Rows are nearest-gene bins, NOT independent signals -- neighbouring genes in one LD block each inherit the block's traits (TTN and TTN-AS1 here are the same 2q31 signal under two symbols). Do not quote the row count as a number of independent loci.",
                          IN_SIG, length(loci)),
  "Mass/wall-thickness block" = sprintf("%s -- still the only mass-primary locus at study-wide, unchanged from the lead-SNP build. At the conventional genome-wide threshold the block had 8 loci (FTO, ACAD10, FGF5, ALDH2, MC4R, PTPN11, ZEB1, LINC00880).",
                          paste(mass_only, collapse = ", ")),
  "Bridging loci" = if (length(bridging))
                      sprintf("%d locus/loci hit both blocks: %s. Placed in the dimension/function block (more hits there) but keeping their mass-block cells colored. MIR4662B is the cross-class finding already flagged as a Word comment on the main text; LINC00964 is NEW at the all-SNP unit and makes the same point.",
                              length(bridging), paste(sapply(bridging, describe_bridge), collapse = "; "))
                    else "none",
  "Cardiac-only loci" = sprintf("%d of %d loci have no study-wide systemic association in the 6 mapped categories: %s.",
                          length(cardiac_only), length(loci), paste(cardiac_only, collapse = ", ")),
  "Empty columns" = sprintf("Echo traits with no study-wide locus: %s. Of these, %s have no genome-wide significant association at all in the source file, so their columns are structurally empty rather than threshold-emptied. Columns are kept to show the absence.",
                          paste(setdiff(CARDIAC, colnames(cardiac_mat)[colSums(cardiac_mat) > 0]), collapse = ", "),
                          if (length(absent)) paste(absent, collapse = ", ") else "none"),
  "Systemic columns" = "Adiposity=Anthropometric, Blood pressure=Vital sign, Metabolic=Metabolism, Liver=Liver, Kidney=Kidney, Blood=Hematology (atlas Category values, 1:1 renamed to match the reference figure's labels).",
  "Cell shade (systemic)" = "0 traits = blank, 1 trait = light orange, >=2 traits = dark orange, all at study-wide significance for that locus.",
  "Method"      = "Input filtered to study-wide significance FIRST, on -log10(P) recomputed from BETA/SE in log space; the file's own study-wide flag column is cross-checked, not trusted. Multi-allelic sites collapsed to one record per SNP-trait pair, matching Table S5/S6.",
  "Panel b"     = sprintf("Unchanged by this revision: LDSC genetic correlations (%s) do not depend on the GWAS significance threshold.", IN_ST7),
  "Script"      = "scripts/echo_partition_rg.R",
  "Generated"   = format(Sys.Date())
))
message("wrote ", OUT_DOC)

# ==== Panel b: echo trait x systemic marker genetic correlation ============
st7 <- normalise_st7(read.csv(IN_ST7, stringsAsFactors = FALSE, check.names = FALSE))

rg_lookup <- function(t1, t2) {
  hit <- st7[(st7$Trait1 == t1 & st7$Trait2 == t2) | (st7$Trait1 == t2 & st7$Trait2 == t1), ]
  if (nrow(hit) == 0) return(c(rg = NA_real_, fdr = NA_real_))
  c(rg = hit$rg[1], fdr = hit$FDR_P[1])
}

rg_mat  <- matrix(NA_real_, nrow = length(PANELB_ROWS), ncol = length(SYS_MARKERS),
                  dimnames = list(PANELB_ROWS, SYS_MARKERS))
fdr_mat <- rg_mat
for (r in PANELB_ROWS) for (c in SYS_MARKERS) {
  v <- rg_lookup(r, c)
  rg_mat[r, c] <- v["rg"]; fdr_mat[r, c] <- v["fdr"]
}
stopifnot(!any(is.na(rg_mat)))

out_b <- data.frame(Trait = PANELB_ROWS, rg_mat, check.names = FALSE)
write.csv(out_b, OUT_B_CSV, row.names = FALSE)
message("wrote ", OUT_B_CSV)
message("panel b: ", sum(fdr_mat < 0.05), " of ", length(fdr_mat), " pairs FDR<0.05; rg range [",
        round(min(rg_mat), 3), ", ", round(max(rg_mat), 3), "]")

# ==== plot ===================================================================
INK        <- "#1a1a1a"
GRID_INK   <- "#e6e6e6"
SEP_INK    <- "#1a1a1a"
BLUE_CELL  <- "#2f5f96"
SYS_LIGHT  <- "#f0b27a"
SYS_DARK   <- "#c8611c"
EMPTY_CELL <- "#f4f4f4"

nL <- length(loci); nCard <- length(CARDIAC); nSys <- length(SYS_COLS)
# Panel a's height now scales with the locus count (33 rows, up from 15).
PANEL_A_H <- max(1150, 34 * nL + 340)
PANEL_B_H <- 1000
png(OUT_PNG, width = 2000, height = PANEL_A_H + PANEL_B_H, res = 150)
layout(matrix(c(1, 2), nrow = 2), heights = c(PANEL_A_H, PANEL_B_H))

## -- panel a --------------------------------------------------------------
par(mar = c(7.5, 12.5, 3.5, 1), family = "sans")

gap <- 1.2
x_card <- seq_len(nCard)
x_sys  <- nCard + gap + seq_len(nSys)
y      <- seq(nL, 1)

plot(NA, xlim = c(0.5, max(x_sys) + 0.5), ylim = c(0.5, nL + 2.2), axes = FALSE, xlab = "", ylab = "")

# background cells
for (j in seq_len(nCard)) rect(x_card[j] - 0.46, y - 0.46, x_card[j] + 0.46, y + 0.46, col = EMPTY_CELL, border = "white")
for (j in seq_len(nSys))  rect(x_sys[j]  - 0.46, y - 0.46, x_sys[j]  + 0.46, y + 0.46, col = EMPTY_CELL, border = "white")

# cardiac fill
for (j in seq_len(nCard)) {
  hit <- cardiac_mat[, j]
  if (any(hit)) rect(x_card[j] - 0.46, y[hit] - 0.46, x_card[j] + 0.46, y[hit] + 0.46, col = BLUE_CELL, border = "white")
}
# systemic fill
for (j in seq_len(nSys)) {
  v <- sys_mat[, j]
  col1 <- ifelse(v == 1, SYS_LIGHT, ifelse(v >= 2, SYS_DARK, NA))
  has <- !is.na(col1)
  if (any(has)) rect(x_sys[j] - 0.46, y[has] - 0.46, x_sys[j] + 0.46, y[has] + 0.46, col = col1[has], border = "white")
}

# vertical separators: mass|dimfn within cardiac block, cardiac|systemic block
abline(v = length(MASS_TRAITS) + 0.5, col = SEP_INK, lwd = 1)
abline(v = (nCard + gap / 2), col = SEP_INK, lwd = 1.6)

# row group separator + row labels (gene names, italic)
grp_break <- sum(grp == "Mass / wall-thickness loci")
if (grp_break > 0 && grp_break < nL) abline(h = y[grp_break] - 0.5, col = SEP_INK, lwd = 1)
row_cex <- if (nL > 26) 0.62 else 0.72
for (i in seq_len(nL)) mtext(loci[i], side = 2, at = y[i], las = 1, line = 0.3, cex = row_cex, font = 4, col = INK)

# group labels (rotated, left of gene names)
if (grp_break > 0) text(par("usr")[1] - 3.6, mean(y[seq_len(grp_break)]), "Mass / wall-thickness loci",
                        srt = 90, cex = 0.78, font = 2, col = INK, xpd = TRUE)
if (grp_break < nL) text(par("usr")[1] - 3.6, mean(y[(grp_break + 1):nL]), "Dimension / function loci",
                         srt = 90, cex = 0.78, font = 2, col = INK, xpd = TRUE)

# column labels (rotated)
for (j in seq_len(nCard)) text(x_card[j], 0.5 - 0.35, CARDIAC[j], srt = 55, adj = 1, cex = 0.78, font = 4, col = INK, xpd = TRUE)
for (j in seq_len(nSys))  text(x_sys[j],  0.5 - 0.35, SYS_COLS[j], srt = 55, adj = 1, cex = 0.78, font = 2, col = INK, xpd = TRUE)

# column group headers
text(mean(x_card), nL + 1.9, "Cardiac traits", cex = 0.95, font = 2, col = INK, xpd = TRUE)
text(mean(x_sys), nL + 1.9, "Systemic traits", cex = 0.95, font = 2, col = SYS_DARK, xpd = TRUE)
segments(min(x_card) - 0.4, nL + 1.55, max(x_card) + 0.4, nL + 1.55, col = INK, xpd = TRUE)
segments(min(x_sys) - 0.4, nL + 1.55, max(x_sys) + 0.4, nL + 1.55, col = INK, xpd = TRUE)

title(main = sprintf("a  Echocardiographic loci partition by trait class (study-wide significant SNPs, %d loci)", nL),
      cex.main = 1.05, col.main = INK, adj = 0, line = 2.1)

## -- panel b --------------------------------------------------------------
par(mar = c(6.2, 12.5, 3.5, 6))

RG_CAP <- 0.8   # symmetric color scale bound
ramp <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
rg_color <- function(v) ramp[pmin(101, pmax(1, round((pmin(pmax(v, -RG_CAP), RG_CAP) + RG_CAP) / (2 * RG_CAP) * 100) + 1))]

nR <- length(PANELB_ROWS); nC <- length(SYS_MARKERS)
xb <- seq_len(nC); yb <- seq(nR, 1)

plot(NA, xlim = c(0.5, nC + 0.5), ylim = c(0.5, nR + 1.3), axes = FALSE, xlab = "", ylab = "")
for (i in seq_len(nR)) for (j in seq_len(nC)) {
  rect(xb[j] - 0.48, yb[i] - 0.48, xb[j] + 0.48, yb[i] + 0.48, col = rg_color(rg_mat[i, j]), border = "white")
  if (fdr_mat[i, j] < 0.05) text(xb[j], yb[i], "*", col = if (abs(rg_mat[i, j]) > RG_CAP * 0.55) "white" else INK,
                                 cex = 1.0, font = 2)
}

b_break <- length(PANELB_MASS)
abline(h = yb[b_break] - 0.5, col = SEP_INK, lwd = 1)
for (i in seq_len(nR)) mtext(PANELB_ROWS[i], side = 2, at = yb[i], las = 1, line = 0.3, cex = 0.78, font = 2, col = INK)
text(par("usr")[1] - 2.3, mean(yb[seq_len(b_break)]), "Mass /\nwall-thickness", srt = 90, cex = 0.78, font = 2, col = INK, xpd = TRUE)
text(par("usr")[1] - 2.3, mean(yb[(b_break + 1):nR]), "Function", srt = 90, cex = 0.78, font = 2, col = INK, xpd = TRUE)

for (j in seq_len(nC)) text(xb[j], 0.5 - 0.42, SYS_MARKERS[j], srt = 55, adj = 1, cex = 0.8, font = 2, col = INK, xpd = TRUE)

title(main = "b  Genetic correlation with systemic markers  (* FDR < 0.05)",
      cex.main = 1.05, col.main = INK, adj = 0, line = 2.1)

# color-scale legend
legend_x <- nC + 1.0
leg_y <- seq(1, nR, length.out = 101)
rect(legend_x, leg_y[-101], legend_x + 0.45, leg_y[-1], col = ramp, border = NA, xpd = TRUE)
rect(legend_x, leg_y[1], legend_x + 0.45, leg_y[101], border = INK, xpd = TRUE)
text(legend_x + 0.7, c(1, (1 + nR) / 2, nR), c(sprintf("-%.1f", RG_CAP), "0.0", sprintf("%.1f", RG_CAP)),
     adj = 0, cex = 0.72, col = INK, xpd = TRUE)
text(legend_x + 0.05, nR + 0.7, expression(r[g]), adj = 0, cex = 0.85, font = 2, col = INK, xpd = TRUE)

dev.off()
message("wrote ", OUT_PNG)
