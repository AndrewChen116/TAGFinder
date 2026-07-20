# TAGFinder Tutorial

This tutorial demonstrates the complete TAGFinder workflow using the example
RNA expression dataset included in this repository.

The example dataset contains:

- 200 samples
- 1,000 genes
- 59 control samples
- 141 cancer samples

The commands below include only required arguments. All optional settings use
the defaults defined by the current TAGFinder scripts. Run the commands from
the TAGFinder repository root.

## Step 0. Prepare the environment

Follow the [installation guideline](https://github.com/AndrewChen116/TAGFinder#installation) to download TAGFinder and setup the environment

> Activate the TAGFinder Conda environment:

```bash
conda activate tagfinder_env
```


Create the result directories automatically by running the modules below. Each
module creates its own output directory if it does not already exist.

## Step 1. Construct phenotype-specific networks

> Module 1 separates samples by phenotype, applies feature quality control, and
constructs correlation, distance, and percentile-ranked feature networks for
the control and cancer groups.

Inputs:

- `example_data/expression.tsv`: samples in rows and genes in columns
- `example_data/metadata.tsv`: sample identifiers and phenotype labels

Run:

```bash
Rscript module1_network_construction.R \
    --feature-table example_data/expression.tsv \
    --metadata example_data/metadata.tsv \
    --phenotype1 control \
    --phenotype2 cancer \
    --outdir results/module1
```

Outputs:

> The default prefix is `analysis`. 

```text
results/module1/
├── analysis_control_correlation_matrix.tsv
├── analysis_control_distance_matrix.tsv
├── analysis_control_topPct_matrix.tsv
├── analysis_control_topPct_numeric.tsv
├── analysis_cancer_correlation_matrix.tsv
├── analysis_cancer_distance_matrix.tsv
├── analysis_cancer_topPct_matrix.tsv
├── analysis_cancer_topPct_numeric.tsv
├── analysis_module1_selected_features.tsv
├── analysis_module1_feature_qc.tsv
├── analysis_module1_matched_samples.tsv
├── analysis_module1_network_summary.tsv
├── analysis_module1_output_manifest.tsv
└── analysis_module2_ripser_input_list.txt
```

The numeric `topPct` matrices contain no row or column names and are passed to
Ripser++ in Module 2. The output manifest retains the labeled matrices required
by later modules.

## Step 2. Perform persistent homology analysis

> Module 2 runs Ripser++ on the two phenotype-specific filtration matrices,
parses the persistence intervals, and summarizes the detected homology
features. By default, homology dimensions up to H1 are calculated, and separate
H0 and H1 topological-index plots are produced with both phenotypes overlaid.

Inputs:

- Module 1 Ripser++ input list
- Ripser++ executable

Run:

```bash
Rscript module2_persistent_homology.R \
    --input-list results/module1/analysis_module2_ripser_input_list.txt \
    --ripser-bin ./ripser-plusplus/ripserplusplus/build/ripser++ \
    --outdir results/module2
```

Outputs:

```text
results/module2/
├── raw_ripser_output/
├── barcode_tables/
├── plots/
│   ├── analysis_module2_H0_topological_index.png
│   ├── analysis_module2_H0_topological_index.pdf
│   ├── analysis_module2_H1_topological_index.png
│   └── analysis_module2_H1_topological_index.pdf
├── analysis_module2_barcode_all.tsv
├── analysis_module2_homology_summary.tsv
├── analysis_module2_plot_manifest.tsv
├── analysis_module2_run_manifest.tsv
└── analysis_module3_input_manifest.tsv
```
<img width="864" height="432" alt="analysis_module2_barcode_plot" src="https://github.com/user-attachments/assets/2f3e71fd-62a5-4fdb-b991-b055ff147b7c" />

Each barcode record contains the homology dimension, birth threshold, death
threshold, lifespan, and infinite-death status of one topological feature.

## Step 3. Determine topological feature composition

> Module 3 integrates the labeled Module 1 network matrices with the Module 2
barcodes to identify the biological features composing each H1 topological
feature. H1 is used by default.

Inputs:

- Module 1 output manifest
- Module 2 manifest prepared for Module 3

Run:

```bash
Rscript module3_feature_composition.R \
    --module1-manifest results/module1/analysis_module1_output_manifest.tsv \
    --module2-manifest results/module2/analysis_module3_input_manifest.tsv \
    --outdir results/module3
```

Main outputs include:

```text
results/module3/
├── composition_tables/
├── analysis_module3_composition_all.tsv
├── analysis_module3_composition_summary.tsv
├── analysis_module3_matched_input_manifest.tsv
└── analysis_module4_input_manifest.tsv
```

The combined composition table links each persistence interval to its component
features and provides the input for participation-score calculation.

## Step 4. Calculate participation scores

> Module 4 calculates feature-level participation scores from the topological
feature compositions. It then compares cancer against control and calculates
the change in participation score for every feature.

Inputs:

- Module 3 manifest prepared for Module 4

Run:

```bash
Rscript module4_participation_score.R \
    --module3-manifest results/module3/analysis_module4_input_manifest.tsv \
    --outdir results/module4 \
    --prefix analysis
```

Outputs:

```text
results/module4/
├── analysis_module4_feature_composition_long.tsv
├── analysis_module4_participation_score_long.tsv
├── analysis_module4_participation_score_wide.tsv
├── analysis_module4_delta_participation_score.tsv
├── analysis_module4_top_abs_delta_participation_score.tsv
├── analysis_module4_score_summary.tsv
└── analysis_module5_input_manifest.tsv
```

For each feature, the delta score is calculated as the cancer participation
score minus the control participation score. A positive value therefore means
increased participation in cancer-associated topological features.

## Step 5. Identify topological-altering genes

> Module 5 tests whether each observed delta participation score is unusual
relative to features with similar network degree. Empirical P-values are
adjusted using the Benjamini-Hochberg procedure to identify topological-altering
genes (TAGs).

Inputs:

- Module 4 manifest prepared for Module 5
- Module 1 output manifest containing the phenotype-specific networks

Run:

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

With the default `topPct` degree matrix, `--degree-threshold 5` defines an edge
as belonging to the top 5% of network connections.

Outputs:

```text
results/module5/
├── plots/
│   ├── analysis_module5_volcano_plot.png
│   └── analysis_module5_volcano_plot.pdf
├── analysis_module5_network_degree_table.tsv
├── analysis_module5_merged_score_degree_table.tsv
├── analysis_module5_degree_matched_null_result.tsv
├── analysis_module5_topological_altering_features.tsv
├── analysis_module5_volcano_data.tsv
├── analysis_module5_summary.tsv
└── analysis_module5_output_manifest.tsv
```

Figure:

<img width="3600" height="1800" alt="analysis_module5_volcano_plot" src="https://github.com/user-attachments/assets/c4f56822-8849-414a-a1a4-00df5f979d50" />

The main final result is:

```text
results/module5/analysis_module5_topological_altering_features.tsv
```

Important result columns include:

| Column | Description |
|---|---|
| `feature_id` | Biological feature identifier, such as a gene symbol |
| `score_tested` | Delta participation score tested in Module 5 |
| `degree_reference` | Network degree used for background matching |
| `n_degree_matched` | Number of degree-matched background features |
| `degree_matched_z` | Standardized delta score relative to the matched background |
| `empirical_p_value` | Degree-matched empirical P-value |
| `q_value` | Benjamini-Hochberg adjusted P-value |
| `direction` | Increased or decreased participation from control to cancer |
