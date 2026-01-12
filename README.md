
## Explainable machine learning reveals drivers of amplifications and deletions across cancer genomes
Fabio Alfieri, Gökçe Senger, Gabriele Oliveto, Manjunatha Kogenaru, Teresa Davoli, Martin Schaefer

For any questions regarding the code and/or data, please contact with Fabio Alfieri (fabio.alfieri@nyulangone.org) or Gabriele Oliveto (gabriele.oliveto@ieo.it).


### PENNE Browser

[**PENNE web browser**](https://schaeferlab.shinyapps.io/PENNE/) — Prediction & Explanation of Non-neutral copy Number Events.

<img src="https://raw.githubusercontent.com/fabio-alfieri/copy-number-annotation/main/PENNE-logo.png" alt="PENNE-logo" width="120"/>


<table>
<tr>
<td>

PENNE web browser available [here](https://schaeferlab.shinyapps.io/PENNE/)
(Prediction & Explanation of Non neutral copy Number Events)

</td>
<td style="width:200px; text-align:right;">

<img src="https://raw.githubusercontent.com/fabio-alfieri/copy-number-annotation/main/PENNE-logo.png"
     alt="PENNE-logo"
     width="200"/>

</td>
</tr>
</table>

### Abstract
Amplifications and deletions of genomic regions are pervasive features of cancer genomes, yet it remains largely unclear which of these focal copy number alterations (CNAs) or chromosome- and arm-level aneuploidies act as drivers of carcinogenesis and which merely reflect underlying chromosomal instability. In this study, we develop an explainable machine learning framework that predicts amplification and deletion frequencies across 11 cancer genomes by integrating genomic-structural features that shape the probability of CNA occurrence with gene-level features indicative of selection. The models achieve high performances across focal, mid-length, arm-level and chromosome-level events, revealing scale-, chromosome- and tumor-dependent selective and mechanistic forces. Local architectural features such as proximity to centromeres, telomeres, and fragile sites mainly drive focal CNAs, whereas mid-length and large-scale CNAs reflect a mixture of structural constraints, dosage sensitivity, and gene-specific selection linked to oncogenes, tumor suppressors, and essential genes. Using SHAP-based interpretability, we generate a genome-wide map that distinguishes regions whose copy-number states are best explained by selective pressure from those arising primarily through structural susceptibility, ultimately providing a web-based annotation browser, named PENNE (Prediction & Explanation of Non-neutral copy Number Events), to investigate the CNA landscape at different scales. Finally, longitudinal single-cell DNA-sequencing of Reversine-induced chromosomal-instability experiments validates the model: early aneuploidies are stochastic, but over time, chromosome arms predicted to confer selective advantage become preferentially retained. Together, our findings establish a framework for interpreting the CNA and aneuploidy landscape of cancer and for systematically uncovering their likely functional drivers.

<img src="https://raw.githubusercontent.com/fabio-alfieri/copy-number-annotation/main/figure-GitHub.png"
     alt="paper-workflow"/>

### System requirements
Analyses were run on Linux:
R version: 
Required R packages (supplementary file with packages link here)


### Reproduce the analyses

#### Clone the repository
```bash
git clone https://github.com/fabio-alfieri/copy-number-annotation.git
```
#### Download supplementary data from Zenodo
(change Zenodo repo link)
```bash
cd path/to/GitHub/copy-number-annotation/
wget -O data.zip https://zenodo.org/record/7079304/files/data.zip?download=1
unzip data.zip
```

#### Reproduce the feature matrix
The R scripts needed to reproduce the feature matrix are stored in
```bash
cd path/to/GitHub/copy-number-annotation/publication-scripts/feature-matrix/compute-features/
```
Scripts are numbered in execution order.

#### Calculate the CNA scores 
The R scripts needed to calculate CNA scores are stored in
```bash
cd path/to/GitHub/copy-number-annotation/publication-scripts/feature-matrix/classify-CNAs/
```
Scripts are numbered in execution order.

#### Run the ML models
The R scripts needed to train the ML models and the subsequent analyses are stored in
```bash
cd path/to/GitHub/copy-number-annotation/publication-scripts/machine-learning/
```
The main script 
