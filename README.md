# scRNAseq-SMC-Has3

This repository contains the full single-cell RNA-seq analysis published in our manuscript:

[Hartmann F, Gorski DJ, _et al._ SMC-Derived Hyaluronan Modulates Vascular SMC Phenotype in Murine Atherosclerosis. Circ Res. 2021 Nov 12;129(11):992-1005.](https://doi.org/10.1161/CIRCRESAHA.120.318479)

Sequencing data can be found at GEO accession [GSM5537130](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM5537130) and [GSM5537131](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM5537131)

## Abstract

<p align="center">
  <img src="/schematic_overview_edit.png" width="1000">
</p>

<p style='text-align: justify;'>
RATIONALE:  Plaque  instability  remains  poorly  understood  and  new  therapeutic  approaches  to  reduce  plaque  rupture  and  subsequent clinical events are of great interest. Recent studies revealed an important role of phenotypic switching of smooth muscle cells (SMC) in controlling plaque stability, including ECM (extracellular matrix) deposition.

<p style='text-align: justify;'>
OBJECTIVE: The aim of this study was to elucidate the role of hyaluronan derived from SMC–hyaluronan synthase 3 (Has3), in phenotypic switching and plaque stability in an animal model of atherosclerosis.

<p style='text-align: justify;'>
METHODS  AND  RESULTS:  A  mouse  line  with  SMC-specific  deletion  of  Has3  and  simultaneous  SMC-lineage  tracing  (eYFP [enhanced yellow fluorescent protein]) on an Apoe−/− background was used. Lineage tracing of SMC with eYFP revealed that SMC-specific deletion of Has3 significantly increased the number of LGALS3+ (galectin-3) transition state SMC and decreased ACTA2+ (alpha-smooth muscle actin) SMC. Notably, SMC-Has3 deletion led to significantly increased collagen deposition  and  maturation  within  the  fibrous  cap  and  the  whole  lesion,  as  evidenced  by  picrosirius  red  staining  and  LC-PolScope  analysis.  Single-cell  RNA  sequencing  of  brachiocephalic  artery  lesions  demonstrated  that  the  loss  of  SMC-Has3 enhanced the transition of SMC to a Lgals3+, ECM-producing phenotype with elevated acute-phase response gene expression. Experiments using cultured murine aortic SMC revealed that blocking CD44 (cluster of differentiation-44), an important hyaluronan binding receptor, recapitulated the enhanced acute-phase response, and synthesis of fibrous ECM.
<p style='text-align: justify;'>
CONCLUSIONS: These studies provide evidence that the deletion of SMC-Has3 results in an ECM-producing transition state SMC phenotype (characterized by LGALS3+ expression), likely via reduced CD44 signaling, resulting in increased collagen formation and maturation, an index consistent with increased plaque stability.

## Analysis
The following analysis was performed with [Seurat](https://satijalab.org/seurat/index.html) v3.1.5, it can be installed with the following command:

```R
remotes::install_version("Seurat", version = "3.1.5")
```

To automatically recreate the published analysis and results, clone this repository and place the count matrices (downloaded from the links above) inside the `data` folder. It should then contain the following:

```
scRNAseq-SMC-Has3/data
  KO_filtered_feature_bc_matrix
  WT_filtered_feature_bc_matrix
  gene_signatures.csv
  published_embeddings.Rdata
```

Running `0-SMC-Has3-complete-analysis.R` will create the necessary results directories and run each analysis step in order. Otherwise, each script can be run individually, starting from `1-SMC-Has3-clustering.R`

<p align="justify">
  <img src="/DimPlot.png" width="24%">
  <img src="/DotPlot_top10_genotype.png" width="75%">
</p>

For consistency, we have incorporated the published UMAP embeddings in the workflow. However, if you would like to process the data on your own, the cells included after pre-processing, their embeddings, annotations and genotypes can be found in the `final_cells.csv` file.
