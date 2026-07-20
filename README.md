# TAGFinder
A bioinformatics framework for feature-level topological attribution in biological networks using Persistent Homology.

---

## Overview

Persistent homology has become an increasingly popular approach for characterizing the higher-order topology of biological networks. Existing persistent homology methods can quantify topological features, such as connected components and cycles, but they cannot determine which biological features are responsible for the observed topological changes.

This software provides an end-to-end framework that identifies feature-level contributions to network topology remodeling by integrating network construction, persistent homology analysis, topological feature composition analysis, participation score calculation, and statistical identification of topological-altering features.

The software is designed for any feature-by-sample matrix, including transcriptomics, DNA methylation, proteomics, metabolomics, and other omics datasets.

<img width="841" height="422" alt="圖片" src="https://github.com/user-attachments/assets/12bc454a-7d96-4ad7-82f2-d86068cb7e65" />


---

## Overview

See [Tutorial.md](Tutorial.md) for the QuickStart workflow.

---

## Software Architecture

| Module | Description |
|---------|-------------|
| 1 | Construct phenotype-specific feature networks |
| 2 | Perform persistent homology analysis |
| 3 | Extract topological feature compositions |
| 4 | Calculate feature participation scores |
| 5 | Identify topological-altering features |

---

## Installation

### Requirements

- Git
- Conda
- CUDA (≥ 10.1)

### 1. Clone repository

```bash
git clone https://github.com/AndrewChen116/TAGFinder.git
cd TAGFinder
```

### 2. Create the Conda environment

The recommended method is to create the environment from environment.yml:
```bash
conda env create -f environment.yml
```

Activate the environment:
```bash
conda activate tagfinder_env
```

### 3. Install Ripser++

Check that the NVIDIA driver and CUDA compiler are available:
```
nvidia-smi
nvcc --version
```

Clone and compile Ripser++ through `install_ripser.sh` :
```
# Make the script executable
chmod +x install_ripser.sh

# Run the automation script
./install_ripser.sh
```

### 4. Validate the TAGFinder installation

Check the command-line interface of each module:
```
Rscript module1_network_construction.R --help
Rscript module2_persistent_homology.R --help
Rscript module3_feature_composition.R --help
Rscript module4_participation_score.R --help
Rscript module5_identify_topological_altering_features.R --help
```

---

## Input

TAGFinder requires two tab-separated input files (tsv):

    1. A feature table containing quantitative measurements.
    2. A metadata table defining the phenotype of each sample.

Sample identifiers must be unique and consistent between the two files.

Example datasets are available in `example_data`


### Feature table

The feature table should contain samples as rows and biological features as columns. The first column must contain sample identifiers, while all remaining columns must contain numeric values.

The biological features may represent genes, CpG sites, proteins, metabolites, or other molecular measurements.

Example (`expression.tsv`, truncated):

```text
sample_id           SFTPC       SCGB1A1    SFTPA1      SFTPA2
TCGA-44-3396-11A    13.9204     10.3068    13.5363     13.9022
TCGA-55-6971-11A    14.9036      1.1403    14.9022     15.5834
TCGA-44-2662-11A    15.0188      3.9638    14.2709     14.3635
TCGA-44-6148-11A    14.2357     12.8331    13.6373     13.8117
TCGA-50-5939-11A    14.4647      7.0348    13.7262     13.9106
```

The example file contains 1,000 genes X 200 samples. Only four features and five samples are shown above.

### Metadata

The metadata table must contain two columns:

- `sample_id`: sample identifiers matching those in the feature table.
- `phenotype`: the phenotype assigned to each sample.

Example (`metadata.tsv`, truncated):

```text
sample_id           phenotype
TCGA-44-3396-11A    control
TCGA-55-6971-11A    control
TCGA-44-2662-11A    control
TCGA-91-6828-01A    cancer
TCGA-91-A4BD-01A    cancer
TCGA-73-4676-01A    cancer
```

The example metadata contains 200 samples: 59 control samples and 141 cancer samples. 




---

## Usage

TAGFinder is a sequential workflow. Run Modules 1–5 in order because each
module creates the standardized inputs or manifest required by the next module.
All commands use the `--argument value` format. Boolean arguments accept
`true`/`false` (case-insensitive).

The examples below use the files in `example_data/`, the default output prefix
`analysis`, and the phenotype comparison `cancer - control`.

### Module 1. Construct phenotype-specific feature networks

> Match the feature table to the metadata, perform feature-level quality
> control, and construct correlation, distance, and percentile-ranked network
> matrices separately for two phenotypes.

**Inputs**

- A TSV/CSV feature table with samples in rows and biological features in columns.
- A TSV/CSV metadata table containing sample identifiers and phenotype labels.

**Usage**

```bash
Rscript module1_network_construction.R \
    --feature-table example_data/expression.tsv \
    --metadata example_data/metadata.tsv \
    --phenotype1 control \
    --phenotype2 cancer \
    --outdir results/module1
```

**Required parameters**

| Parameter | Value | Description |
|---|---|---|
| `--feature-table` | File path | Sample-by-feature measurement table in TSV or CSV format. Rows represent samples and columns represent biological features. |
| `--metadata` | File path | Sample metadata table containing the sample identifier and phenotype columns. |
| `--phenotype1` | String | Baseline phenotype label. It must exactly match a value in the metadata phenotype column. |
| `--phenotype2` | String | Comparison phenotype label. It must differ from `--phenotype1` and exactly match a metadata value. |
| `--outdir` | Directory path | Directory in which Module 1 outputs are written. It is created automatically if absent. |

**Optional parameters**

| Parameter | Accepted values | Default | Description |
|---|---|---|---|
| `--prefix` | String | `analysis` | Prefix added to all Module 1 output filenames. |
| `--sample-col` | Column name | `sample_id` | Sample identifier column in the metadata table. |
| `--phenotype-col` | Column name | `phenotype` | Phenotype column in the metadata table. |
| `--feature-id-col` | Column name | First non-numeric column, otherwise row names | Column in the feature table containing sample identifiers. Despite the historical argument name, this identifies rows/samples rather than biological features. |
| `--cor-method` | `pearson`, `spearman` | `pearson` | Correlation method used to construct each phenotype-specific feature network. |
| `--distance-transform` | `one_minus_abs_cor`, `one_minus_cor` | `one_minus_abs_cor` | Converts correlation to distance. `one_minus_abs_cor` uses `1 - abs(r)`; `one_minus_cor` uses `(1 - r) / 2`. |
| `--feature-filter-mode` | `common_valid`, `none` | `common_valid` | `common_valid` retains features passing missingness, nonzero, and variance criteria in both phenotypes; `none` disables this filtering. |
| `--min-samples-per-phenotype` | Integer ≥ 1 | `3` | Minimum number of matched samples required in each phenotype. |
| `--min-nonmissing-fraction` | Numeric in `[0,1]` | `0.8` | Minimum fraction of non-missing observations required for a feature in each phenotype. |
| `--min-nonzero-fraction` | Numeric in `[0,1]` | `0` | Minimum fraction of nonzero observations required for a feature in each phenotype. |
| `--impute-missing` | `none`, `median` | `none` | Missing-value strategy. With `none`, correlations use pairwise-complete observations; `median` imputes each feature within a phenotype before correlation. |
| `--round-digits` | Non-negative integer | `6` | Number of digits retained in the exported network matrices. |
| `--write-correlation` | `true`, `false` | `true` | Whether to write the phenotype-specific correlation matrices. |

**Main outputs**

- Labeled correlation, distance, and `topPct` matrices for both phenotypes.
- Header-free numeric `topPct` matrices for Ripser++.
- Feature QC, matched-sample, selected-feature, and network-summary tables.
- `<prefix>_module1_output_manifest.tsv` for Modules 3 and 5.
- `<prefix>_module2_ripser_input_list.txt` for Module 2.

---

### Module 2. Perform persistent homology analysis

> Run Ripser++ on the phenotype-specific filtration matrices, parse persistence
> intervals, summarize homology features, and generate overlaid H0/H1
> topological-index plots.

**Inputs**

- The Ripser++ matrix-path list generated by Module 1.
- A compiled Ripser++ executable.

**Usage**

```bash
Rscript module2_persistent_homology.R \
    --input-list results/module1/analysis_module2_ripser_input_list.txt \
    --ripser-bin ./ripser-plusplus/ripserplusplus/build/ripser++ \
    --outdir results/module2
```

**Required parameters**

| Parameter | Value | Description |
|---|---|---|
| `--input-list` | File path | Text file generated by Module 1. Each non-empty, non-comment line must contain one numeric filtration-matrix path. |
| `--ripser-bin` | Executable path | Path to the compiled Ripser++ executable. This argument remains required in parse-only mode. |
| `--outdir` | Directory path | Directory in which Module 2 outputs are written. |

**Optional parameters**

| Parameter | Accepted values | Default | Description |
|---|---|---|---|
| `--prefix` | String | `analysis` | Prefix added to all Module 2 output filenames. |
| `--format` | Ripser++ format string | `distance` | Matrix format passed to the Ripser++ `--format` argument. |
| `--max-dim` | Integer ≥ 0 | `1` | Maximum homology dimension calculated by Ripser++. The default calculates H0 and H1; the plotting function currently visualizes H0 and H1 only. |
| `--threshold` | Numeric ≥ 0 | `100` | Maximum filtration threshold passed to Ripser++. |
| `--ratio` | Numeric > 0 | `1` | Approximation ratio passed to Ripser++. |
| `--sparse` | `true`, `false` | `true` | Whether to add the Ripser++ `--sparse` flag. |
| `--gpu-id` | CUDA device ID | Not set | Sets `CUDA_VISIBLE_DEVICES` for the Ripser++ process. |
| `--extra-args` | Quoted string | Empty | Additional Ripser++ arguments appended after the standard arguments. |
| `--validate-matrix` | `true`, `false` | `true` | Validate that every input is a square, header-free numeric matrix before execution. |
| `--run-ripser` | `true`, `false` | `true` | Set to `false` to parse existing raw outputs without rerunning Ripser++. Expected raw files must already exist at the standard output paths. |
| `--continue-on-error` | `true`, `false` | `true` | Continue processing the remaining matrices if one Ripser++ run fails. |
| `--make-plots` | `true`, `false` | `true` | Generate overlaid H0/H1 topological-index plots. H0 is indexed by death and H1 is indexed by birth. |
| `--plot-format` | `png`, `pdf`, `both` | `both` | Output format for the topological-index plots. |
| `--plot-labels` | Comma-separated labels | Matrix IDs | Legend labels in input-list order. The number of labels must equal the number of matrices. |
| `--plot-colors` | Comma-separated R colors | `#56B4E9,#E64B4B` | Group colors used in the overlaid plots. Each group is drawn with 50% opacity. |
| `--plot-width` | Numeric > 0 | `6` | Figure width in inches. |
| `--plot-height` | Numeric > 0 | `6` | Figure height in inches. |
| `--plot-dpi` | Integer > 0 | `300` | PNG resolution in dots per inch. This does not affect PDF output. |

**Main outputs**

- Raw Ripser++ text outputs and per-matrix barcode tables.
- `<prefix>_module2_barcode_all.tsv` and homology summary.
- H0/H1 topological-index figures and a plot manifest.
- `<prefix>_module2_run_manifest.tsv` with commands, exit codes, and errors.
- `<prefix>_module3_input_manifest.tsv` for Module 3.

---

### Module 3. Extract topological feature compositions

> Match persistence intervals to the labeled filtration matrices and identify
> which biological features compose each target topological feature.

**Inputs**

- The Module 1 output manifest containing labeled filtration matrices.
- The Module 2 manifest containing barcode tables.

**Usage**

```bash
Rscript module3_feature_composition.R \
    --module1-manifest results/module1/analysis_module1_output_manifest.tsv \
    --module2-manifest results/module2/analysis_module3_input_manifest.tsv \
    --outdir results/module3
```

**Required parameters**

| Parameter | Value | Description |
|---|---|---|
| `--module1-manifest` | File path | `<prefix>_module1_output_manifest.tsv` generated by Module 1. Provides phenotype labels and labeled `topPct` matrices. |
| `--module2-manifest` | File path | `<prefix>_module3_input_manifest.tsv` generated by Module 2. Provides matrix IDs and barcode-table paths. |
| `--outdir` | Directory path | Directory in which Module 3 outputs are written. |

**Optional parameters**

| Parameter | Accepted values | Default | Description |
|---|---|---|---|
| `--prefix` | String | `analysis` | Prefix added to all Module 3 output filenames. |
| `--target-dim` | Integer ≥ 0 | `1` | Homology dimension for composition extraction. The implementation and biological interpretation are primarily designed for H1. |
| `--mode` | `birth_death`, `birth_only` | `birth_death` | `birth_death` extracts both birth and death compositions; `birth_only` skips death-composition extraction. |
| `--birth-choice` | `largest`, `first`, `smallest` | `largest` | Strategy for resolving multiple candidate boundary pairs at the same birth value, based on candidate composition size or encounter order. |
| `--value-digits` | Integer ≥ 1 | `10` | Significant digits used to create numeric matching keys between barcode thresholds and matrix values. |
| `--max-boundary-pairs` | Integer ≥ 1 | `20000` | Safety cap on candidate matrix boundary pairs retained for each filtration value. |
| `--keep-failed` | `true`, `false` | `true` | Keep intervals whose composition could not be resolved, with diagnostic status fields. |
| `--verbose` | `true`, `false` | `true` | Print per-matrix progress messages. |

**Main outputs**

- Per-phenotype composition tables in `composition_tables/`.
- `<prefix>_module3_composition_all.tsv` with birth, death, and total composition.
- Composition summary and matched-input manifest.
- `<prefix>_module4_input_manifest.tsv` for Module 4.

---

### Module 4. Calculate participation scores

> Attribute topological features to their component biological features,
> calculate phenotype-specific participation scores, and compute
> `phenotype2 - phenotype1` delta scores.

**Inputs**

Provide either the Module 3 combined composition table or its Module 4 input
manifest. Using the manifest is recommended for the sequential workflow.

**Usage**

```bash
Rscript module4_participation_score.R \
    --module3-manifest results/module3/analysis_module4_input_manifest.tsv \
    --outdir results/module4 \
    --prefix analysis
```

**Required parameters**

| Parameter | Value | Description |
|---|---|---|
| `--outdir` | Directory path | Directory in which Module 4 outputs are written. |
| `--prefix` | String | Prefix added to all Module 4 output filenames. |

**Required input (choose one)**

| Parameter | Value | Description |
|---|---|---|
| `--module3-manifest` | File path | Module 3 manifest containing the `composition_all` path. Provide this or `--composition`. |
| `--composition` | File path | Explicit Module 3 `<prefix>_module3_composition_all.tsv` path. Provide this or `--module3-manifest`; if both are supplied, this explicit path takes precedence. |

**Optional parameters**

| Parameter | Accepted values | Default | Description |
|---|---|---|---|
| `--phenotype1` | Phenotype label | First phenotype in the composition table | Baseline phenotype used in delta calculation. |
| `--phenotype2` | Phenotype label | Second phenotype in the composition table | Comparison phenotype. Delta scores are calculated as `phenotype2 - phenotype1`. |
| `--composition-col` | `total_composition`, `birth_composition`, `death_composition` | `total_composition` | Composition column used to attribute topological features to biological features. |
| `--feature-list` | TSV file path | Not set | Optional complete feature list. Features absent from all compositions can be added with zero participation scores. |
| `--feature-col` | Column name | `feature_id` | Feature identifier column in `--feature-list`. |
| `--min-lifespan` | Numeric ≥ 0 | `0` | Minimum persistence lifespan retained for participation-score calculation. |
| `--include-infinite-death` | `true`, `false` | `false` | Include intervals with infinite death instead of excluding them. |
| `--infinite-lifespan-value` | Numeric | `0` | Finite lifespan assigned to infinite-death intervals when they are included. |
| `--separator` | String | `;` | Delimiter used between feature IDs in composition columns. |
| `--write-zero-scores` | `true`, `false` | `true` | Add zero-score rows for known features absent from the selected compositions. |
| `--top-n` | Integer ≥ 1 | `100` | Number of features with the largest absolute delta scores written to the top-delta table. |

**Main outputs**

- Long-form feature composition and participation-score tables.
- Wide phenotype-specific participation-score table.
- `<prefix>_module4_delta_participation_score.tsv`.
- Top absolute-delta table and score summary.
- `<prefix>_module5_input_manifest.tsv` for Module 5.

---

### Module 5. Identify topological-altering features

> Test whether each delta participation score is unusual relative to features
> with similar network degree, adjust empirical P-values, identify
> topological-altering features, and generate a volcano plot.

**Inputs**

- The Module 1 output manifest containing phenotype-specific network matrices.
- Either the Module 4 delta-score table or the Module 5 input manifest generated
  by Module 4.

**Usage**

```bash
Rscript module5_identify_topological_altering_features.R \
    --module4-manifest results/module4/analysis_module5_input_manifest.tsv \
    --module1-manifest results/module1/analysis_module1_output_manifest.tsv \
    --outdir results/module5 \
    --prefix analysis \
    --phenotype1 control \
    --phenotype2 cancer \
    --degree-threshold 5
```

**Required parameters**

| Parameter | Value | Description |
|---|---|---|
| `--module1-manifest` | File path | `<prefix>_module1_output_manifest.tsv` generated by Module 1. Provides the phenotype-specific matrices used for degree calculation. |
| `--outdir` | Directory path | Directory in which Module 5 outputs are written. |
| `--prefix` | String | Prefix added to all Module 5 output filenames. |
| `--phenotype1` | Phenotype label | Baseline phenotype. Must match a phenotype in the Module 1 manifest and the Module 4 delta table. |
| `--phenotype2` | Phenotype label | Comparison phenotype. Must differ from `--phenotype1`. Positive delta scores represent increased participation in this phenotype. |
| `--degree-threshold` | Numeric | Network edge threshold. For `topPct`, values must be in `[0,100]` and an edge satisfies `0 < topPct <= threshold`; for `distance`, an edge satisfies `0 < distance <= threshold`. |

**Required input (choose one)**

| Parameter | Value | Description |
|---|---|---|
| `--module4-manifest` | File path | Module 4 manifest containing the delta participation-score path. Provide this or `--delta-score`. |
| `--delta-score` | File path | Explicit Module 4 delta participation-score table. Provide this or `--module4-manifest`; if both are supplied, this explicit path takes precedence. |

**Optional parameters**

| Parameter | Accepted values | Default | Description |
|---|---|---|---|
| `--degree-matrix-type` | `topPct`, `distance` | `topPct` | Module 1 matrix type used to define network edges and calculate degree. |
| `--score-col` | Numeric column name | `delta_participation_score_lifespan_sum` | Delta participation-score column tested against the degree-matched background. |
| `--degree-reference` | `avg`, `phenotype1`, `phenotype2` | `avg` | Degree used for matching. `avg` is the mean degree across both phenotypes. |
| `--min-matched` | Integer ≥ 1 | `20` | Minimum number of degree-matched background features required for an empirical test. |
| `--max-window` | Numeric ≥ 0 | `50` | Largest permitted absolute degree difference when expanding the matching window. |
| `--alternative` | `two.sided`, `greater`, `less` | `two.sided` | Alternative hypothesis for the empirical test. The two-sided test uses deviation from the matched-score median. |
| `--q-cutoff` | Numeric in `[0,1]` | `0.05` | Benjamini-Hochberg q-value cutoff for formal TAF calls. |
| `--p-cutoff` | Numeric in `[0,1]`, `NA` | `NA` | Optional additional raw empirical P-value cutoff for TAF calls. `NA` disables this filter. |
| `--min-abs-delta` | Numeric ≥ 0 | `0` | Minimum absolute tested delta score required for a formal TAF call. |
| `--exclude-zero-degree` | `true`, `false` | `false` | Exclude degree-zero features from empirical testing. |
| `--pseudocount` | Numeric ≥ 0 | `1` | Pseudocount added to the numerator and denominator of empirical P-value calculations. |
| `--make-volcano` | `true`, `false` | `true` | Generate the Module 5 volcano plot. |
| `--volcano-format` | `png`, `pdf`, `both` | `both` | Volcano-plot output format. |
| `--volcano-p-cutoff` | Numeric in `(0,1]` | `0.05` | Raw empirical P-value threshold used for volcano-plot coloring and the horizontal reference line. This is independent of the q-value cutoff for formal TAF calls. |
| `--volcano-delta-cutoff` | Numeric ≥ 0 | Value of `--min-abs-delta` | Absolute delta-score threshold used for volcano-plot coloring and vertical reference lines. |
| `--volcano-width` | Numeric > 0 | `12` | PNG/PDF width in inches. The wider default accommodates the right-side legend. |
| `--volcano-height` | Numeric > 0 | `6` | PNG/PDF height in inches. |
| `--volcano-dpi` | Integer > 0 | `300` | PNG resolution in dots per inch. This does not affect PDF output. |
| `--volcano-point-size` | Numeric > 0 | `0.65` | Point-size multiplier used in the volcano plot. |

**Main outputs**

- Network degree and merged score-degree tables.
- Degree-matched null-test results with empirical P-values and BH q-values.
- `<prefix>_module5_topological_altering_features.tsv` containing formal TAF calls.
- Volcano plot data plus PNG/PDF figures.
- Module 5 summary, run information, and output manifest.

---

## Output

Final outputs include

- Participation scores
- Delta participation scores
- Degree-matched statistics
- Topological-altering genes (TAGs)

---

## Contact

Kuan-Lin Chen

Institute of Plant and Microbial Biology, Academia Sinica

Email: kuanlin@gate.sinica.edu.tw
