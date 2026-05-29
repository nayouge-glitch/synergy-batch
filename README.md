# SynergyFinder batch runner

Compute **ZIP, Bliss, HSA, and Loewe** synergy scores for *all* your dose-response
matrices in one run, using the official `synergyfinder` R package — the same engine
as the SynergyFinder web app. Because this runs on GitHub's servers (full internet
access), it installs `drc` and produces numbers that match the website.

## How to use

1. Create a new GitHub repository and upload these files, keeping the structure:

   ```
   run_synergy.R
   .github/workflows/synergy.yml
   input/        <- put your matrix .xlsx files here
   output/       <- results appear here (as a downloadable artifact)
   ```

2. Drop your matrix files into `input/`. Use the SynergyFinder matrix layout
   (the same format you already generate), with **% inhibition** values:

   ```
   Drug1:    | Lapatinib
   Drug2:    | Trastuzumab
   ConcUnit: | nM
            |  0    0.1   1    10   100
   0        |  ...
   0.1      |  ...
   ...
   (blank row, then the next block)
   ```

   You can stack many blocks per file and many files — all are processed together.

3. Push to GitHub. The workflow runs automatically (or trigger it from the
   **Actions** tab → *Run workflow*).

4. When it finishes, open the run and download the **`synergy-results`** artifact.
   It contains:
   - `synergy_summary.xlsx` — one row per block with mean ZIP/Bliss/HSA/Loewe.
   - `synergy_scores_full.xlsx` — per-dose-pair scores (for landscapes).
   - `synergy_result.rds` — full R object for further analysis/plots.

## Matching the website exactly

- Readout is set to **inhibition** in `run_synergy.R`.
- Baseline correction is set to `"all"` (variable `CORRECT_BASELINE`). If your
  website runs used **Correction OFF**, change it to `"non"`. Options: `non`, `part`, `all`.
- Curve fitting is LL.4 (the package default), same as the web app.

## Notes

- First run installs packages and may take ~15–20 min; later runs are cached.
- Loewe is undefined for pairs where one drug alone has ~0% effect (e.g. an
  antibody with no single-agent activity); the package returns NA there.
