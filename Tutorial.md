# Tutorial

This tutorial demonstrates the complete workflow of TAGFinder using the example RNA expression dataset.

The example dataset contains

- 200 samples
- 200 genes
- phenotype1 (100 samples)
- phenotype2 (100 samples)

Running this tutorial will reproduce the complete analysis pipeline and identify Topological Altering Genes (TAGs).

---

# Step 0. Prepare files

Project structure

```

TAGFinder/
│
├── module1_network_construction.R
├── module2_persistent_homology.R
├── module3_feature_composition.R
├── module4_participation_score.R
├── module5_identify_topological_altering_features.R
│
├── example/
│   ├── expression.tsv
│   └── metadata.tsv
│
└── results/

```

---

# Step 1. Construct feature-feature networks

Input

```

example/expression.tsv

example/metadata.tsv

```

Run

```bash
Rscript module1_network_construction.R \
    --feature-table example/expression.tsv \
    --metadata example/metadata.tsv \
    --sample-col sample_id \
    --phenotype-col phenotype \
    --phenotype1 phenotype1 \
    --phenotype2 phenotype2 \
    --outdir results/module1 \
    --prefix example
```

Expected outputs

```

results/module1/

phenotype1/
correlation.tsv
distance.tsv
topPct.tsv
topPct.numeric.tsv

phenotype2/
correlation.tsv
distance.tsv
topPct.tsv
topPct.numeric.tsv

network_summary.tsv

module1_manifest.tsv

module2_ripser_input_list.txt

```

---

# Step 2. Persistent Homology Analysis

Input

```

results/module1/module2_ripser_input_list.txt

```

Run

```bash
Rscript module2_persistent_homology.R \
    --input-list results/module1/module2_ripser_input_list.txt \
    --ripser-bin /path/to/ripser++ \
    --outdir results/module2 \
    --prefix example \
    --max-dim 1
```

Expected outputs

```

results/module2/

barcode/

homology_summary.tsv

module2_run_manifest.tsv

module3_input_manifest.tsv

```

---

# Step 3. Topological Feature Composition Analysis

Run

```bash
Rscript module3_feature_composition.R \
    --module1-manifest results/module1/module1_manifest.tsv \
    --module2-manifest results/module2/module3_input_manifest.tsv \
    --outdir results/module3 \
    --prefix example \
    --target-dim 1
```

Expected outputs

```

results/module3/

composition_all.tsv

composition_summary.tsv

module4_input_manifest.tsv

```

---

# Step 4. Participation Score Calculation

Run

```bash
Rscript module4_participation_score.R \
    --module3-manifest results/module3/module4_input_manifest.tsv \
    --outdir results/module4 \
    --prefix example
```

Expected outputs

```

results/module4/

participation_score_long.tsv

participation_score_wide.tsv

delta_participation_score.tsv

module5_input_manifest.tsv

```

---

# Step 5. Identify Topological Altering Genes

Run

```bash
Rscript module5_identify_topological_altering_features.R \
    --module4-manifest results/module4/module5_input_manifest.tsv \
    --module1-manifest results/module1/module1_manifest.tsv \
    --outdir results/module5 \
    --prefix example \
    --phenotype1 phenotype1 \
    --phenotype2 phenotype2 \
    --degree-matrix-type topPct \
    --degree-threshold 5
```

Expected outputs

```

results/module5/

network_degree_table.tsv

degree_matched_null_result.tsv

topological_altering_features.tsv

summary.tsv

```

---

# Final Results

The final output file is

```

results/module5/

example_module5_topological_altering_features.tsv

```

Each row corresponds to one candidate Topological Altering Gene (TAG).

The output includes

| Column | Description |
|----------|-------------|
| feature_id | Gene name |
| participation_score_phenotype1 | Participation score in phenotype1 |
| participation_score_phenotype2 | Participation score in phenotype2 |
| delta_participation_score | Difference between phenotypes |
| degree_phenotype1 | Network degree in phenotype1 |
| degree_phenotype2 | Network degree in phenotype2 |
| empirical_pvalue | Degree-matched empirical P-value |
| qvalue | Benjamini-Hochberg adjusted P-value |
| significant | Significant TAG |

---

# Workflow Summary

```

expression.tsv
metadata.tsv

↓

Module 1
Network Construction

↓

Module 2
Persistent Homology

↓

Module 3
Composition Analysis

↓

Module 4
Participation Score

↓

Module 5
TAG Identification

↓

Topological Altering Genes

```

---

# Running the complete pipeline

Alternatively, execute all modules sequentially.

```bash
Rscript module1_network_construction.R ...
Rscript module2_persistent_homology.R ...
Rscript module3_feature_composition.R ...
Rscript module4_participation_score.R ...
Rscript module5_identify_topological_altering_features.R ...
```
