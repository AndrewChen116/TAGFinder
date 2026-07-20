# TAGFinder
A bioinformatics framework for feature-level topological attribution in biological networks using Persistent Homology.

---

## Overview

Persistent homology has become an increasingly popular approach for characterizing the higher-order topology of biological networks. Existing persistent homology methods can quantify topological features, such as connected components and cycles, but they cannot determine which biological features are responsible for the observed topological changes.

This software provides an end-to-end framework that identifies feature-level contributions to network topology remodeling by integrating network construction, persistent homology analysis, topological feature composition analysis, participation score calculation, and statistical identification of topological-altering features.

The software is designed for any feature-by-sample matrix, including transcriptomics, DNA methylation, proteomics, metabolomics, and other omics datasets.

<img width="841" height="422" alt="圖片" src="https://github.com/user-attachments/assets/12bc454a-7d96-4ad7-82f2-d86068cb7e65" />


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

- R (≥ 4.3)
- Git
- CUDA >=10.1

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

### Feature table

Rows represent biological samples.

Columns represent biological features.

```
Sample    Feature1    Feature2    Feature3
S1
S2
S3
...
```

---

### Metadata

```
Sample    Phenotype

S1        phenotype1

S2        phenotype1

S3        phenotype2
```

---

## Usage

### Module 1

```bash
Rscript module1_network_construction.R \
    --feature-table expression.tsv \
    --metadata metadata.tsv
```

---

## Module 2

```bash
Rscript module2_persistent_homology.R \
    --module1-manifest module1_manifest.tsv
```

---

### Module 3

```bash
Rscript module3_feature_composition.R \
    --module2-manifest module2_manifest.tsv
```

---

## Module 4

```bash
Rscript module4_participation_score.R \
    --module3-manifest module3_manifest.tsv
```

---

### Module 5

```bash
Rscript module5_identify_topological_altering_features.R \
    --module4-manifest module4_manifest.tsv
```

---

## Output

```
results/

module1/

module2/

module3/

module4/

module5/

summary/

logs/
```

Final outputs include

- Participation scores
- Delta participation scores
- Degree-matched statistics
- Topological-altering features (TAFs)

---

## Example Dataset

Example datasets are available in

```
example/

feature_table.tsv

metadata.tsv
```

Example command

```bash
bash run_example.sh
```


---

## Contact

Kuan-Lin Chen

Institute of Plant and Microbial Biology, Academia Sinica

Email: kuanlin@gate.sinica.edu.tw
