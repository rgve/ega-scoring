# A Universal, library‐aware Quality Score for NGS Files

A fast, **data.table**-based pipeline to compute a **continuous QC score** per sequencing file, **tailored to library strategy**.  
The workflow: **load → merge → summarize (robust) → score → export → (optional) plot**.

> This repository contains an R script that ingests a metrics table and a mapping from EGAF IDs to library strategies, computes robust strategy-specific thresholds, and assigns per-file QC scores without categorical labels.

## Features

- **Library-strategy aware** thresholds and scoring
- **Robust stats** per (filetype × metric × strategy): `median`, `IQR`, `P5`, `P95`, `MAD`
- **Hard-fail** bump if a metric value is outside `[P5, P95]`
- **Continuous score** (sum of metric penalties), no categorical buckets
- Optional **boxplots** and **histograms** for top strategies
- Scales to large datasets using **data.table**


## Requirements
```
- **R ≥ 4.1**
- Packages:
  - [`data.table`](https://rdatatable.io/)
  - [`ggplot2`](https://ggplot2.tidyverse.org/) *(optional; for plots)*
  - [`scales`](https://scales.r-lib.org/) *(optional; for nicer plot axes)*
```
Install packages:

```r
install.packages(c("data.table", "ggplot2", "scales"))
```

 *(If you won’t plot, `ggplot2` and `scales` are not required—set `plot_enable <- FALSE`.)*

## Inputs (Data Schemas)

### 1) `parsed_all_120k.csv` — **metrics table**
- **Format:** CSV with **header**
- **Row grain:** one row per `(egaf_id, filetype, flag_key)` observation
- **Required columns:**
  - `egaf_id` *(string)* — e.g., `EGAF00003178367`
  - `filetype` *(string)* — expected values include `fastq`, `bamcram`, `vcf`, `__error__`
  - `flag_key` *(string)* — metric name. **This pipeline uses numeric metrics:**
    - For `fastq`: `gc_content`, `duplicate_reads`, `quality_reads`
    - For `bamcram`: `gc_content`, `duplicate_reads`, `mapq`, `unaligned`
  - `value` *(numeric or coercible string)*
- **Example:**
  ```csv
  egaf_id,filetype,flag_key,value
  EGAF00000001,fastq,gc_content,42.1
  EGAF00000001,fastq,duplicate_reads,9.8
  EGAF00000001,fastq,quality_reads,83.5
  EGAF00000002,bamcram,gc_content,44.2
  EGAF00000002,bamcram,duplicate_reads,12.7
  EGAF00000002,bamcram,mapq,96.0
  EGAF00000002,bamcram,unaligned,5.1
  ```

### 2) `file_all_lib_results.csv` — **EGAF → library strategy mapping**
- **Format:** CSV **without header** (4 columns)
- **Row grain:** one row per EGAF
- **Columns (in order):**
  1. `egaf_id` *(string)*
  2. `run_stable_id` *(string)* — not used by the pipeline
  3. `experiment_stable_id` *(string)* — not used by the pipeline
  4. `library_strategy` *(string)* — e.g., `WGS`, `WGA`, `RNA-Seq`, `ATAC-seq`, `ChIP-Seq`
- **Example (no header):**
  ```csv
  EGAF00000001,ERR000001,EGAX000001,WGS
  EGAF00000002,ERR000002,EGAX000002,RNA-Seq
  EGAF00000003,ERR000003,EGAX000003,ATAC-seq
  ```

> Rows missing a strategy (unmapped `egaf_id`) will be dropped when restricting to **top-K** strategies per filetype.


## How it works

1. **Load** metrics and strategy mapping; merge by `egaf_id`.
2. **Filter** to numeric metrics needed by filetype; coerce `value` → numeric.
3. **Select top-K** `library_strategy` per `filetype` to avoid sparse tails.
4. **Summarize** value distributions per `(filetype, flag_key, library_strategy)`:
   - `N`, `median`, `IQR`, `P5`, `P95`, `MAD`.
   - Floor very small `MAD` to the **5th percentile of positive MADs** to stabilize scores.
5. **Score** each observation with strategy-aware rules (examples):
   - `gc_content`: `abs(value - median) / MAD`
   - `duplicate_reads`: `pmax(0, (value - median)/MAD)`
   - `quality_reads` (FASTQ): `pmax(0, (80 - value)/10)`
   - `mapq` (BAM/CRAM): `pmax(0, (80 - (100 - value))/10)`
   - **Clamp** each metric score to `[0, 5]`
   - **Hard fail**: if `value < P5` or `value > P95`, bump that metric score to **≥ 2**
6. **Aggregate** per file: `quality_score = sum(metric_score)`.

No categorical labels are produced—**consumers set their own cut-offs**.


## Configuration

Inside the script you’ll find these knobs:

```r
top_k_per_filetype <- 20                          # # strategies per filetype to keep
numeric_keys_fastq <- c("duplicate_reads","gc_content","quality_reads")
numeric_keys_bam   <- c("duplicate_reads","gc_content","mapq","unaligned")
plot_enable        <- TRUE                        # set FALSE to skip plotting
seed_sampling      <- 42                          # reproducible plots
```

**Paths to edit:**
```r
input_file   <- "/path/to/parsed_all_120k.csv"
lib_map_file <- "/path/to/file_all_lib_results.csv"
out_dir      <- "/path/to/qc_outputs"
```

## Outputs

### `qc_outputs/summary_stats_topk.csv`
Robust thresholds per `(filetype, flag_key, library_strategy)`:
- Columns: `filetype, flag_key, library_strategy, N, median, IQR, P5, P95, MAD`

### `qc_outputs/final_scores.csv`
Continuous score per file:
- Columns: `egaf_id, filetype, library_strategy, quality_score`

### Plots (optional)
- **Boxplots**: `qc_outputs/plots_box_top20/box_<filetype>_<flag_key>.png`
- **Histograms**: `qc_outputs/plots_hist_top20/hist_<filetype>_<flag_key>.png`


## Validation checks (recommended)

Add these before scoring if desired:

```r
stopifnot(all(c("egaf_id","filetype","flag_key","value") %in% names(dt)))
stopifnot(all(c("egaf_id","library_strategy") %in% names(lib_map)))
if (anyNA(data_full$library_strategy)) message("Note: some EGAFs lack a mapped library_strategy.")
```


## Performance tips

- Keep `top_k_per_filetype` moderate (e.g., 10–30) to maintain robust thresholds.
- Use `data.table` keys for joins (already done).
- Avoid plotting on extremely large datasets or turn off with `plot_enable <- FALSE`.

## Troubleshooting

- **All zero/NA scores?** Check that `flag_key`s match the expected names and `value` is numeric.
- **Exploding scores?** Ensure `MAD` flooring is active (step 4b).
- **Missing strategies?** Verify `egaf_id` consistency between the two CSVs.
