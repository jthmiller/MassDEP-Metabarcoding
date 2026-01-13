# MassDEP-Metabarcoding

## RBCL Metabarcoding of MA streams 

## Processing Sequecnces
1. code/pipeline.sh
- remove polyg tails with 
- match and trim rbcl primer sequences with qiime2::cutadapt
- denoise with qiime2::dada2

2. 

## TITAN Analysis 
See [Smucker et al.](https://pubmed.ncbi.nlm.nih.gov/32602216/), [Pilgrim et al.](https://pubmed.ncbi.nlm.nih.gov/36213613/)
Threshold indicator taxa analysis [TITAN2](https://cran.r-project.org/web/packages/TITAN2/vignettes/titan2-intro.pdf)
- "OTUs that increased or decreased along TP and TN gradients along with nutrient concentrations at which assemblages had substantial changes in the occurrences and relative abundances of OTUs."
- "threshold indicator taxa analysis (TITAN), boosted regression tree analyses, and gradient forest analysis: nonparametric statistical analyses characterized diatom OTU and assemblage relationships with TP and TN concentrations"
-  "TITAN incorporates relative abundance and occurrence data to quantify multi‐taxa assemblage change, but only handles one predictor variable at a time"
- "NMDS, we used 50 runs of real data and 1,000 randomizations (PC‐ORD v. 5, MjM Software, Gleneden Beach, Oregon, USA)."
## Plots





## Notes

https://academic.oup.com/ismecommun/article/5/1/ycaf171/8265810?login=false
(assignTaxonomy() function), with a minimum bootstrap value of 60% and the ready to use Diat.barcode reference library v.12 for metabarcoding analyses (doi.org/10.57745/XWJJGI) which is an adaptation of the original library
- n = 156, stop codon were detected by converting ASV nucleotide sequences into amino acid sequences using Emboss v6.6.0), bringing the total number of ASVs across the 1103 samples from 12 404 to 12 248. From 12 248 ASVs, 8650 remained after filtering out sequences that could not be assigned at the order level, and further reduced to 3302 ASVs across 1073 samples by removing rare ASVs (fewer than 10 sequences in total) and samples with fewer than 500 reads

Biostrings:: codons(x)


cubar::seq_to_codons(yeast_cds_qc[["YDR320W-B"]])