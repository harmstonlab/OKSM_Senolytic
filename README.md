# OKSM_Senolytic

## Versions
STAR: 2.7.1a
RSEM: v1.3.1
AME: 5.1.0
FASTQC: v0.11.8
R:
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0  (64-bit)
Running under: macOS Big Sur 11.3.1
Versions of R packages used can be found in the Session Info section at the bottom of each .Rmd file.

## Data
All raw and aligned data for this analysis can be obtained from GEO via the study accession number XXXXXX. The processed data can also be found in this repository under the folder `data`.
Annotation file used: Drosophila_melanogaster.BDGP6.22.97.chr.gtf (Ensembl v97)
Reference fasta file used: dm6.fa (Ensembl v97)

## Data labels
RNA-seq label|Description
---|---
Control| Guts from control flies.
O (OKSM)|Guts from flies in which OKSM is over-expressed in the intestinal stem cells.
S (Senolytic)|Guts from flies with the senolytic construct expressed in the intestinal stem cells.
SO (Senolytic_OKSM)|Guts from flies in which the senolytic construct is expressed and OKSM is over-expressed in the intestinal stem cells.
TdTom| Guts from control flies.
OKSM| Guts from flies in which OKSM is over-expressed in the ubiquitous expression model.
Sen (Senolytic)|Guts from flies with the senolytic construct expressed in the ubiquitous expression model.
Sen_OKSM (Senolytic_OKSM)|Guts from flies in which the senolytic construct is expressed and OKSM is over-expressed in the ubiquitous expression model.

## Analysis
**! Ensure that the annotation file is in a folder called `annotation` before running the analysis**

The analysis was carried out to compare gene expression between 4 conditions -- control/TdTomato (Control), OKSM (Ab), Senolytic (Sen), and OKSM and Senolytic together (Sen_OKSM) -- across two datasets, using the isc.Rmd and ubiquitous.Rmd for the analysis. 
