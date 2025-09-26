# ────────────────────────────────────────────────────────────────────────────────
# Continuous scoring workflow: load → merge → summarize → score
# ────────────────────────────────────────────────────────────────────────────────

# Libraries ---------------------------------------------------------------------
# data.table is used throughout for speed and clear keyed joins.
# ggplot2/scales are only needed if plot_enable == TRUE.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)  # optional plots
  library(scales)   # pretty axes for plots
})

# Inputs (edit paths) -----------------------------------------------------------
# parsed_all_120k.csv: red-flag metrics (wide variety of flags; we keep numeric ones only)
# file_all_lib_results.csv: EGAF → library_strategy mapping (no header; 4 columns)
input_file   <- "/Users/raul/Desktop/parsed_all_120k.csv"
lib_map_file <- "/Users/raul/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/ega.nosync/bioteam/infer_hd_report/file_all_lib_results.csv"

# Outputs -----------------------------------------------------------------------
# All CSVs and plots are written under this folder.
out_dir         <- "/Users/raul/Desktop/qc_outputs"
out_summary_csv <- file.path(out_dir, "summary_stats_topk.csv")  # thresholds by (type×metric×strategy)
out_scores_csv  <- file.path(out_dir, "final_scores.csv")        # per-file continuous score
out_box_dir     <- file.path(out_dir, "plots_box_top20")         # optional boxplots
out_hist_dir    <- file.path(out_dir, "plots_hist_top20")        # optional histograms
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters --------------------------------------------------------------------
# Top-K strategies per filetype keeps thresholds representative (avoids sparse tails).
top_k_per_filetype <- 20

# Metrics we consider numeric and score-able by filetype:
numeric_keys_fastq <- c("duplicate_reads","gc_content","quality_reads")
numeric_keys_bam   <- c("duplicate_reads","gc_content","mapq","unaligned")

# Toggle plotting. Set to FALSE for a compute-only run.
plot_enable        <- TRUE

# Seed for reproducibility in any sampling step inside the plotting loop.
seed_sampling      <- 42

# 0) Load data ------------------------------------------------------------------
# Read the red-flag table (expects header)
dt <- fread(input_file)

# Read mapping (no header); name columns explicitly
lib_map_raw <- fread(
  lib_map_file, header = FALSE,
  col.names = c("egaf_id","run_stable_id","experiment_stable_id","library_strategy")
)

# 1) Merge library strategy -----------------------------------------------------
# Keep only the columns needed downstream for the join and reporting
lib_map <- lib_map_raw[, .(egaf_id, library_strategy)]

# Key the tables by egaf_id for a fast merge
setkey(dt, egaf_id); setkey(lib_map, egaf_id)

# Left-join the strategy into the metrics (rows with unknown egaf_id get NA strategy)
data_full <- merge(dt, lib_map, by = "egaf_id", all.x = TRUE)

# 2) Keep numeric metrics & coerce value ----------------------------------------
# Define the universe of numeric metrics we will score
numeric_keys <- unique(c(numeric_keys_fastq, numeric_keys_bam))

# Filter: drop VCF and internal error rows; keep only the chosen numeric flag_keys.
# Coerce 'value' to numeric (invalid strings become NA).
data_num <- data_full[
  !filetype %in% c("vcf","__error__") & flag_key %in% numeric_keys
][
  , value := suppressWarnings(as.numeric(value))
][]

# 3) Select top-K library strategies per filetype -------------------------------
# Count rows by (filetype, library_strategy) to rank strategies by prevalence.
counts <- data_num[, .N, by = .(filetype, library_strategy)][order(filetype, -N)]

# For each filetype, keep the top-K strategies (ties broken deterministically by order)
topk   <- counts[, head(.SD, top_k_per_filetype), by = filetype][, .(filetype, library_strategy)]

# Join this filter onto the numeric data to keep only those strategies
setkey(topk, filetype, library_strategy); setkey(data_num, filetype, library_strategy)
data_topk <- data_num[topk, nomatch = 0]

# 4) Robust per-strategy summaries (ONLY what’s used downstream) ----------------
# We compute robust summaries per (filetype × flag_key × library_strategy):
#  - N:   number of finite values
#  - median, IQR, P5, P95, MAD: used in scoring and hard-fail thresholds
summary_stats_topk <- data_topk[, .(
  N      = sum(!is.na(value)),
  median = stats::median(value, na.rm = TRUE),
  IQR    = IQR(value, na.rm = TRUE),
  P5     = stats::quantile(value, 0.05, na.rm = TRUE, names = FALSE),
  P95    = stats::quantile(value, 0.95, na.rm = TRUE, names = FALSE),
  MAD    = mad(value, na.rm = TRUE)
), by = .(filetype, flag_key, library_strategy)][order(filetype, flag_key, -N)]

# 4b) Stabilize very small MAD values ------------------------------------------
# Very small MADs can explode scores; we floor MAD to the 5th percentile of positive MADs.
mad_vec <- summary_stats_topk$MAD
min_mad <- if (any(mad_vec > 0, na.rm = TRUE)) {
  stats::quantile(mad_vec[mad_vec > 0], 0.05, na.rm = TRUE, names = FALSE)
} else { 1.0 }
summary_stats_topk[, MAD := fifelse(MAD < min_mad, min_mad, MAD)]

# 5) Strategy-aware scoring -----------------------------------------------------
# Prepare the thresholds table and join onto each observation in data_topk.
thr <- summary_stats_topk[, .(filetype, flag_key, library_strategy, median, IQR, P5, P95, MAD)]
setkey(thr, filetype, flag_key, library_strategy)
setkey(data_topk, filetype, flag_key, library_strategy)
dt_merge <- thr[data_topk, nomatch = 0]

# Compute per-metric raw scores (metric_score_raw) + hard-fail flags in one pass.
# Then clamp to [0,5] and apply hard-fail bump to ≥2.
dt_merge[
  , `:=`(
    metric_score_raw = fcase(
      # FASTQ metrics
      filetype == "fastq"   & flag_key == "gc_content",      abs((value - median) / MAD),
      filetype == "fastq"   & flag_key == "duplicate_reads", pmax(0, (value - median) / MAD),
      filetype == "fastq"   & flag_key == "quality_reads",   pmax(0, (80 - value) / 10),

      # BAM/CRAM metrics
      filetype == "bamcram" & flag_key == "gc_content",      abs((value - median) / MAD),
      filetype == "bamcram" & flag_key == "duplicate_reads", pmax(0, (value - median) / MAD),
      filetype == "bamcram" & flag_key == "unaligned",       pmax(0, (value - median) / MAD),
      filetype == "bamcram" & flag_key == "mapq",            pmax(0, (80 - (100 - value))/10),

      # Fallback
      default = as.numeric(NA)
    ),
    # Hard-fail if outside the robust central band
    hard = (value < P5 | value > P95)
  )
][
  # Clamp each metric score to [0,5]
  , metric_score := pmin(pmax(metric_score_raw, 0), 5)
][
  # Ensure hard-fail metrics contribute at least 2 points of penalty
  hard == TRUE, metric_score := pmax(metric_score, 2)
][]

# 6) Aggregate to file-level continuous score -----------------------------------
# Sum per-metric scores across metrics available for each file.
# (No categorical labels; downstream consumers can choose cutoffs appropriate to their use.)
final_scores <- dt_merge[
  , .(quality_score = sum(metric_score, na.rm = TRUE)),
  by = .(egaf_id, filetype, library_strategy)
][order(filetype, -quality_score)]

# 7) Export ---------------------------------------------------------------------
# Save thresholds (for auditing/tuning) and per-file scores (for consumers).
fwrite(summary_stats_topk, out_summary_csv)
fwrite(final_scores,      out_scores_csv)

# 8) Optional plots -------------------------------------------------------------
# For each (filetype, flag_key), draw:
#  - boxplots per top-20 library_strategy (faceted 2x10, free y-scale)
#  - histograms of value distributions per strategy
if (isTRUE(plot_enable)) {
  dir.create(out_box_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(out_hist_dir, showWarnings = FALSE, recursive = TRUE)

  combos <- unique(data_topk[, .(filetype, flag_key)])
  set.seed(seed_sampling)

  for (i in seq_len(nrow(combos))) {
    ft <- combos$filetype[i]; fk <- combos$flag_key[i]

    # Subset and sanity check
    sub <- data_topk[filetype == ft & flag_key == fk & is.finite(value)]
    if (nrow(sub) == 0) next

    # Identify the top-20 strategies for this combo (by count)
    counts_c  <- sub[, .N, by = library_strategy][order(-N)]
    top20     <- counts_c[1:min(20, .N), library_strategy]

    # Filter to those strategies and set facet order
    sub20     <- sub[library_strategy %in% top20]
    sub20[, library_strategy := factor(library_strategy, levels = top20)]
    dt_counts <- counts_c[library_strategy %in% top20][, library_strategy := factor(library_strategy, levels = top20)]

    # Boxplots (distribution summaries)
    p_box <- ggplot() +
      geom_boxplot(
        data = sub20,
        aes(x = 1, y = value, fill = library_strategy),
        width = 0.7, outlier.shape = 16, outlier.size = 0.6, outlier.alpha = 0.2, alpha = 0.6
      ) +
      geom_text(data = dt_counts, aes(x = 1, label = paste0("n=", N)),
                y = -Inf, vjust = -0.5, size = 3) +
      facet_wrap(~ library_strategy, nrow = 2, ncol = 10, scales = "free_y") +
      scale_y_continuous(labels = comma) +
      labs(title = sprintf("%s — %s", ft, fk),
           subtitle = "Top 20 library strategies (by count)",
           x = NULL, y = "Value") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            strip.text = element_text(size = 8, angle = 45, hjust = 1),
            legend.position = "none", panel.spacing = unit(0.5, "lines"))
    ggsave(file.path(out_box_dir, sprintf("box_%s_%s.png", ft, fk)),
           p_box, width = 16, height = 6, dpi = 300)

    # Histograms (shape of value distributions)
    p_hist <- ggplot(sub20, aes(x = value, fill = library_strategy)) +
      geom_histogram(bins = 30, color = "white", alpha = 0.8) +
      facet_wrap(~ library_strategy, nrow = 2, ncol = 10, scales = "free_y") +
      scale_x_continuous(labels = comma) +
      scale_y_continuous(labels = comma) +
      labs(title = sprintf("%s — %s", ft, fk),
           subtitle = "Value distribution for top 20 library strategies",
           x = "Value", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none",
            strip.text = element_text(size = 8, angle = 45, hjust = 1),
            panel.spacing = unit(0.5, "lines"))
    ggsave(file.path(out_hist_dir, sprintf("hist_%s_%s.png", ft, fk)),
           p_hist, width = 16, height = 6, dpi = 300)
  }
}

# 9) Quick console summary ------------------------------------------------------
# Lightweight descriptive stats of the resulting per-file scores, by filetype.
cat("\n— Summary per filetype (quality_score):\n")
print(final_scores[, .(
  n      = .N,
  mean   = mean(quality_score),
  median = median(quality_score),
  sd     = sd(quality_score),
  IQR    = IQR(quality_score)
), by = filetype][order(filetype)])

cat("\nSaved:\n  • thresholds → ", out_summary_csv,
    "\n  • scores     → ", out_scores_csv, "\n")


# Violin plots: per-metric score distributions, grouped by filetype -------------
# Requires: dt_merge with columns egaf_id, filetype, flag_key, metric_score

library(data.table)
library(ggplot2)

# Keep finite scores; retain strategy in case you want to facet later by it.
df_scores <- as.data.table(dt_merge)[is.finite(metric_score),
                                     .(egaf_id, filetype, flag_key, library_strategy, metric_score)
]

# Enforce a consistent metric ordering on the x-axis.
metric_order <- c("gc_content","duplicate_reads","quality_reads","mapq","unaligned")
df_scores[, flag_key := factor(flag_key, levels = intersect(metric_order, unique(flag_key)))]

# Violin plot with median dots; y-axis limited to the score clamp [0,5].
p_violin <- ggplot(df_scores, aes(x = flag_key, y = metric_score)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = median, geom = "point", size = 0.8) +
  facet_wrap(~ filetype, scales = "free_x", nrow = 2) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold")
  ) +
  labs(
    title = "Per-metric QC score distributions by filetype",
    x = "Metric (flag_key)",
    y = "Per-metric score (0–5)"
  )

print(p_violin)
