# =====================================================================
# run_synergy.R
# Batch synergy scoring with the official SynergyFinder R package.
# Computes ZIP, Bliss, HSA, Loewe for EVERY block in EVERY input file,
# in a single run. Output goes to output/.
#
# Input: any .xlsx in input/ using the SynergyFinder matrix layout
#        (the exampleMatrix-style blocks):
#
#   Drug1:    | Lapatinib
#   Drug2:    | Trastuzumab
#   ConcUnit: | nM
#            |  0    0.1   1    10   100      <- Drug2 concentrations
#   0        |  ...                           <- Drug1 conc, then responses
#   0.1      |  ...
#   ...
#   (blank row, then the next block)
#
# Responses must be % INHIBITION (100 - viability).
# =====================================================================

suppressMessages({
  library(synergyfinder)
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
})

input_dir  <- "input"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

# ---- which baseline correction to match the website -----------------
# The web app "Correction ON" corresponds to correcting all values using
# the fitted single-agent baseline. If your site numbers were produced
# with Correction OFF, set this to "non".
CORRECT_BASELINE <- "non"   # one of: "non", "part", "all"
# ---------------------------------------------------------------------

files <- list.files(input_dir, pattern = "\\.xlsx?$", full.names = TRUE)
if (length(files) == 0) stop("No .xlsx files found in input/")

parse_blocks <- function(path) {
  raw <- suppressMessages(read_excel(path, sheet = 1, col_names = FALSE))
  m <- as.matrix(as.data.frame(raw, stringsAsFactors = FALSE))
  nr <- nrow(m)
  blocks <- list(); i <- 1
  while (i <= nr) {
    a <- trimws(as.character(m[i, 1]))
    if (!is.na(a) && a %in% c("Drug1:", "Drug1")) {
      drug1 <- trimws(as.character(m[i, 2]))
      drug2 <- trimws(as.character(m[i + 1, 2]))
      unit  <- trimws(as.character(m[i + 2, 2]))
      hdr <- suppressWarnings(as.numeric(m[i + 3, -1]))
      hdr <- hdr[!is.na(hdr)]
      ncol2 <- length(hdr)
      r <- i + 4; conc1 <- c(); vals <- list()
      while (r <= nr) {
        c1 <- suppressWarnings(as.numeric(m[r, 1]))
        if (is.na(c1)) break
        conc1 <- c(conc1, c1)
        vals[[length(vals) + 1]] <- suppressWarnings(as.numeric(m[r, 2:(1 + ncol2)]))
        r <- r + 1
      }
      blocks[[length(blocks) + 1]] <- list(
        drug1 = drug1, drug2 = drug2, unit = unit,
        conc1 = conc1, conc2 = hdr, mat = do.call(rbind, vals))
      i <- r
    } else i <- i + 1
  }
  blocks
}

# ---- build one long-format data frame across all files/blocks -------
long <- list(); bid <- 0; meta <- list()
for (f in files) {
  for (b in parse_blocks(f)) {
    bid <- bid + 1
    meta[[bid]] <- data.frame(block_id = bid, source = basename(f),
                              drug1 = b$drug1, drug2 = b$drug2,
                              stringsAsFactors = FALSE)
    for (ii in seq_along(b$conc1)) for (jj in seq_along(b$conc2)) {
      long[[length(long) + 1]] <- data.frame(
        block_id = bid, drug1 = b$drug1, drug2 = b$drug2,
        conc1 = b$conc1[ii], conc2 = b$conc2[jj],
        response = b$mat[ii, jj],
        conc_unit1 = b$unit, conc_unit2 = b$unit,
        stringsAsFactors = FALSE)
    }
  }
}
df   <- do.call(rbind, long)
meta <- do.call(rbind, meta)
cat("Parsed", bid, "blocks from", length(files), "file(s).\n")

# ---- run all four models in one call --------------------------------
data <- ReshapeData(df, data_type = "inhibition")
res  <- CalculateSynergy(
  data,
  method = c("ZIP", "HSA", "Bliss", "Loewe"),
  correct_baseline = CORRECT_BASELINE
)

# ---- summary table: mean synergy per model, per block ---------------
dp <- res$drug_pairs
syn_cols <- grep("_synergy$", colnames(dp), value = TRUE)

if (length(syn_cols) > 0) {
  summary_tbl <- dp[, c("block_id", "drug1", "drug2", syn_cols), drop = FALSE]
} else {
  # fallback: average per-dose-pair scores over combination cells only
  ss <- res$synergy_scores
  model_cols <- grep("_synergy$", colnames(ss), value = TRUE)
  summary_tbl <- ss %>%
    filter(conc1 > 0, conc2 > 0) %>%
    group_by(block_id) %>%
    summarise(across(all_of(model_cols), ~mean(.x, na.rm = TRUE)),
              .groups = "drop")
}
summary_tbl <- left_join(meta, summary_tbl, by = "block_id")

write.xlsx(summary_tbl, file.path(output_dir, "synergy_summary.xlsx"))
write.xlsx(res$synergy_scores, file.path(output_dir, "synergy_scores_full.xlsx"))
saveRDS(res, file.path(output_dir, "synergy_result.rds"))

cat("\n==== SYNERGY SUMMARY (all models) ====\n")
print(summary_tbl, row.names = FALSE)
cat("\nWrote: output/synergy_summary.xlsx, synergy_scores_full.xlsx, synergy_result.rds\n")
