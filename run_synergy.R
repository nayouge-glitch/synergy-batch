# =====================================================================
# run_synergy.R  (v4)
# All four models (ZIP, Bliss, HSA, Loewe) for every block, one run.
# Baseline correction = "part" to match the SynergyFinder web app.
# NO global jitter (v3's jitter perturbed ZIP scores). Data is passed
# through unchanged; each block is wrapped in tryCatch so a single
# failure does not stop the rest.
# Responses must be % INHIBITION.
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

# Set to the SAME option you used on the website.
CORRECT_BASELINE <- "part"   # "non", "part", or "all"

files <- list.files(input_dir, pattern = "\\.xlsx?$", full.names = TRUE)
if (length(files) == 0) stop("No .xlsx files found in input/")

parse_blocks <- function(path) {
  raw <- suppressMessages(read_excel(path, sheet = 1, col_names = FALSE))
  m <- as.matrix(as.data.frame(raw, stringsAsFactors = FALSE))
  nr <- nrow(m); blocks <- list(); i <- 1
  while (i <= nr) {
    a <- trimws(as.character(m[i, 1]))
    if (!is.na(a) && a %in% c("Drug1:", "Drug1")) {
      drug1 <- trimws(as.character(m[i, 2]))
      drug2 <- trimws(as.character(m[i + 1, 2]))
      unit  <- trimws(as.character(m[i + 2, 2]))
      hdr <- suppressWarnings(as.numeric(m[i + 3, -1])); hdr <- hdr[!is.na(hdr)]
      ncol2 <- length(hdr); r <- i + 4; conc1 <- c(); vals <- list()
      while (r <= nr) {
        c1 <- suppressWarnings(as.numeric(m[r, 1]))
        if (is.na(c1)) break
        conc1 <- c(conc1, c1)
        vals[[length(vals) + 1]] <- suppressWarnings(as.numeric(m[r, 2:(1 + ncol2)]))
        r <- r + 1
      }
      blocks[[length(blocks) + 1]] <- list(drug1 = drug1, drug2 = drug2,
        unit = unit, conc1 = conc1, conc2 = hdr, mat = do.call(rbind, vals))
      i <- r
    } else i <- i + 1
  }
  blocks
}

all_blocks <- list(); bid <- 0
for (f in files) for (b in parse_blocks(f)) {
  bid <- bid + 1; b$block_id <- bid; b$source <- basename(f)
  all_blocks[[bid]] <- b
}
cat("Parsed", bid, "blocks from", length(files), "file(s).\n")

run_one <- function(b) {
  M <- b$mat
  M[is.na(M)] <- 0          # only fill NA; do NOT jitter real values
  long <- expand.grid(ii = seq_along(b$conc1), jj = seq_along(b$conc2))
  df <- data.frame(
    block_id = 1, drug1 = b$drug1, drug2 = b$drug2,
    conc1 = b$conc1[long$ii], conc2 = b$conc2[long$jj],
    response = M[cbind(long$ii, long$jj)],
    conc_unit1 = b$unit, conc_unit2 = b$unit, stringsAsFactors = FALSE)

  data <- ReshapeData(df, data_type = "inhibition")
  res <- CalculateSynergy(data, method = c("ZIP","HSA","Bliss","Loewe"),
                          correct_baseline = CORRECT_BASELINE)
  dp <- res$drug_pairs
  syn_cols <- grep("_synergy$", colnames(dp), value = TRUE)
  if (length(syn_cols) > 0) {
    out <- dp[1, syn_cols, drop = FALSE]
  } else {
    ss <- res$synergy_scores; mc <- grep("_synergy$", colnames(ss), value = TRUE)
    out <- as.data.frame(lapply(ss[ss$conc1>0 & ss$conc2>0, mc, drop=FALSE],
                                mean, na.rm = TRUE))
  }
  cbind(data.frame(block_id=b$block_id, source=b$source,
                   drug1=b$drug1, drug2=b$drug2, stringsAsFactors=FALSE), out)
}

rows <- list()
for (b in all_blocks) {
  cat("Block", b$block_id, ":", b$drug1, "+", b$drug2, "... ")
  r <- tryCatch(run_one(b), error = function(e) {
    cat("FAILED:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) { cat("ok\n"); rows[[length(rows)+1]] <- r }
}
if (length(rows) == 0) stop("All blocks failed.")
summary_tbl <- dplyr::bind_rows(rows)
write.xlsx(summary_tbl, file.path(output_dir, "synergy_summary.xlsx"))
cat("\n==== SYNERGY SUMMARY (all models, correction =", CORRECT_BASELINE, ") ====\n")
print(summary_tbl, row.names = FALSE)
