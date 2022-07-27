# OKSM_Senolytic

## Versions
* **STAR:** 2.7.1a
* **RSEM:** v1.3.1
* **AME:** 5.1.0
* **FASTQC:** v0.11.8
* **R:** 
```
* R version 4.1.2 (2021-11-01)
* Platform: x86_64-apple-darwin17.0  (64-bit)
* Running under: macOS Big Sur 11.3.1
* Versions of R packages used can be found in the Session Info section at the bottom of each .Rmd file.
```
## Data
All raw and aligned data for this analysis can be obtained from GEO via the study accession number [GSE201338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201338). The processed data can also be found in this repository under the folder `data`.

Annotation file used: `Drosophila_melanogaster.BDGP6.22.97.chr.gtf` (Ensembl v97)

Reference fasta file used: `dm6.fa` (Ensembl v97)

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

The analysis was carried out to compare gene expression between 4 conditions -- control/TdTomato (Control), OKSM (Ab), Senolytic (Sen), and OKSM and Senolytic together (Sen_OKSM) -- across two datasets, using the `isc.Rmd` and `ubiquitous.Rmd` R Markdown files for the analysis. 

## Construct analysis
To look at how much of the OKSM construct is present in the relevant conditions, we first obtain the sequence of the OKSIM construct from [here](https://www.addgene.org/24603/sequences/). Next, the `.gtf` file was modified to contain 3 extra entries:
```
oksm	AddedGenes	exon	608	12417	.	+	.	gene_id "oksm"; gene_name "oksm"; gene_biotype "construct"; transcript_id "oksm";
oksm	AddedGenes	exon	608	1927	.	+	.	gene_id"myc"; gene_name "myc"; gene_biotype "construct"; transcript_id "myc";
oksm	AddedGenes	exon	8736	12383	.	+	.	gene_id "oks"; gene_name "oks"; gene_biotype "construct"; transcript_id "oks";
```
The locations for these are based on the information obtained from [here](https://www.addgene.org/browse/sequence/172960/).

The necessary files -- `oksm.fa` and the modified `Drosophila_melanogaster.BDGP6.22.97.chr.gtf` -- can be found on the server at `/home/ellora/projects/oksm_nov2021/construct/`

The new reference genome is then generated using STAR and the command:

```
STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /home/ellora/projects/oksm_nov2021/with_constructs/ensembl97 \
--genomeFastaFiles dm6.fa oksm.fa \
--sjdbGTFfile Drosophila_melanogaster.BDGP6.22.97.chr.gtf \
--sjdbOverhang 100 \
--genomeSAindexNbases 13
```

After successfully generating the new genome, the reference is prepared for `rsem`:
```
rsem-prepare-reference --gtf Drosophila_melanogaster.BDGP6.22.97.chr.gtf -p 10 dm6.fa,oksm.fa oksm
```

The script (`oksm_construct.sh`) used to then align the samples to these new references can be found on the server at `/home/ellora/projects/oksm_nov2021/construct`. 