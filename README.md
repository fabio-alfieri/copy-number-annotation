
# Explainable machine learning reveals drivers of amplifications and deletions across cancer genomes
*Fabio Alfieri, Gökçe Senger, Gabriele Oliveto, Manjunatha Kogenaru, Teresa Davoli, Martin Schaefer*

<img align="right"
     src="https://raw.githubusercontent.com/fabio-alfieri/copy-number-annotation/main/PENNE-logo.png"
     width="250"
     alt="PENNE web browser">
This manuscript also presents [PENNE (Prediction & Explanation of Non-neutral copy Number Events)](https://schaeferlab.shinyapps.io/PENNE/), a publicly accessible web platform that supports interactive, pan-cancer exploration of curated copy-number annotation profiles.

For any questions regarding the code and/or data, please contact with Fabio Alfieri (fabio.alfieri@nyulangone.org) or Gabriele Oliveto (gabriele.oliveto@ieo.it).

## Abstract

Amplifications and deletions of genomic regions are pervasive features of cancer genomes, yet it remains largely unclear which of these focal copy number alterations (CNAs) or chromosome- and arm-level aneuploidies act as drivers of carcinogenesis and which merely reflect underlying chromosomal instability. In this study, we develop an explainable machine learning framework that predicts amplification and deletion frequencies across 11 cancer genomes by integrating genomic-structural features that shape the probability of CNA occurrence with gene-level features indicative of selection. The models achieve high performances across focal, mid-length, arm-level and chromosome-level events, revealing scale-, chromosome- and tumor-dependent selective and mechanistic forces. Local architectural features such as proximity to centromeres, telomeres, and fragile sites mainly drive focal CNAs, whereas mid-length and large-scale CNAs reflect a mixture of structural constraints, dosage sensitivity, and gene-specific selection linked to oncogenes, tumor suppressors, and essential genes. Using SHAP-based interpretability, we generate a genome-wide map that distinguishes regions whose copy-number states are best explained by selective pressure from those arising primarily through structural susceptibility, ultimately providing a web-based annotation browser, named PENNE (Prediction & Explanation of Non-neutral copy Number Events), to investigate the CNA landscape at different scales. Finally, longitudinal single-cell DNA-sequencing of Reversine-induced chromosomal-instability experiments validates the model: early aneuploidies are stochastic, but over time, chromosome arms predicted to confer selective advantage become preferentially retained. Together, our findings establish a framework for interpreting the CNA and aneuploidy landscape of cancer and for systematically uncovering their likely functional drivers.

<img src="https://raw.githubusercontent.com/fabio-alfieri/copy-number-annotation/main/figure-GitHub.png"
     alt="paper-workflow"/>


## Reproduce the Analyses 

### Clone GitHub Repository and Download the Data
Clone the GitHub repository on your local machine.
```bash
git clone https://github.com/fabio-alfieri/copy-number-annotation.git
```

**(change Zenodo link with final after publishing repo)**

Download [here](https://zenodo.org/api/records/17737479/draft/files/data.zip/content) (Zenodo) the data.zip folder (required), move the data.zip in the cloned GitHub folder and unzip it.
```bash
cd path/to/GitHub/copy-number-annotation/
wget -O data.zip https://zenodo.org/api/records/17737479/draft/files/data.zip/content
unzip data.zip
```

### Train the XGBoost Models
The R scripts needed to train the ML models and the subsequent analyses are located at `path/to/GitHub/copy-number-annotation/publication-scripts/machine-learning/`
Please run `main_workflow.R` ad follows:
```bash
cd path/to/GitHub/copy-number-annotation/publication-scripts/machine-learning/
Rscript main_workflow.R
```

#### Compute the SHAP-based CNA Annotation
```bash
cd path/to/GitHub/copy-number-annotation/publication-scripts/machine-learning/
Rscript compute_sum_abs_SHAPs.R
```
Annotation matrices and plots are stored at `data/annotation/merged_res_annot_*` and `data/plots/`, respectively.
Annotation matrices are also available on Zenodo [here](https://zenodo.org/api/records/17737479/draft/files/annotation-matrix.zip/content) or:
```bash
cd path/to/GitHub/copy-number-annotation/
wget -O annotation-matrix.zip https://zenodo.org/api/records/17737479/draft/files/annotation-matrix.zip/content
unzip annotation-matrix.zip
```

### Supplementary 
#### Reproduce the Feature Matrix
Pre-built feature matrix are available in .tsv on Zenodo [here](https://zenodo.org/api/records/17737479/draft/files/feature-matrix.zip/content), or as follows:
```bash
cd path/to/GitHub/copy-number-annotation/
wget -O feature-matrix.zip https://zenodo.org/api/records/17737479/draft/files/feature-matrix.zip/content
unzip feature_matrix.zip
```

The R scripts needed to reproduce the feature matrix are located at
```bash
path/to/GitHub/copy-number-annotation/publication-scripts/feature-matrix/compute-features/
```
Scripts are numbered in execution order.





## System requirements
Analyses were run on Linux:
R version: 
Required R packages (supplementary file with packages link here): R-packages.yaml

