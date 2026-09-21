
R
library('phyloseq')
library('qiime2R')
feats <- read_qza('MEALGAE_MERGED_table.qza')$data
seq <- read_qza('MEALGAE_MERGED_rep-seqs.qza')$data
tax <- read_qza('MEALGAE-RBCLNX0923_02122024_Rsyst_vsearch_taxonomy.qza')$data

tax_all <- read_qza('MEALGAE-RBCLNX0923_diat_barcode_v10_263bp-sklearn-hybrid-classification_all.qza')$data
tax_95 <- read_qza('MEALGAE-RBCLNX0923_diat_barcode_v10_263bp-sklearn-hybrid-classification_10_95.qza')$data
tax_90 <- read_qza('MEALGAE-RBCLNX0923_diat_barcode_v10_263bp-sklearn-hybrid-classification_10_90.qza')$data



tax <- read_qza('MEALGAE-RBCLNX0923_02122024_
rownames(tax) <- tax$Feature.ID


sort(feats[,'ME659b071822'])



tax['724a68f3ac4b16c4f976e153247a12d8',]

tax['6f685a09b41c1363ad21bf4ffb972178',]
as.character(seq['6f685a09b41c1363ad21bf4ffb972178'])

  tophit_vsearch_taxonomy.qza
80f40ddc1bb21d5b6d563ed83c4095a6 84675ba9b5147acea952fcdcac7b2b6c 
                           47143                            56916 
0980b245a3d8944c078c2c86c8e7d7a9 f81899aeb313723ddec251c7506e02d8 
                           74221                            76553 
724a68f3ac4b16c4f976e153247a12d8 6f685a09b41c1363ad21bf4ffb972178 
                           91167                           213448 



                                                k__Chromista;p__Bacillariophyta;c__Bacillariophyceae 
                                                                                                   7 
            k__Chromista;p__Bacillariophyta;c__Bacillariophyceae;o__Achnanthales;f__Achnanthidiaceae 
                                                                                                   1 
k__Chromista;p__Bacillariophyta;c__Bacillariophyceae;o__Bacillariales;f__Bacillariaceae;g__Nitzschia 
                                                                                                   2 
                                 k__Chromista;p__Bacillariophyta;c__Bacillariophyceae;o__Naviculales 
                                                                                                   1 
                                 k__Chromista;p__Bacillariophyta;c__Mediophyceae;o__Thalassiosirales 
                                                                                                   1 
                                                                                          Unassigned 
                                                                                                   8 


                                       Eukaryota;Chromista;Chromobiota;Bacillariophyta;Bacillariophyceae 
                                                                                                       8 
         Eukaryota;Chromista;Chromobiota;Bacillariophyta;Bacillariophyceae;Achnanthales;Achnanthidiaceae 
                                                                                                       1 
Eukaryota;Chromista;Chromobiota;Bacillariophyta;Bacillariophyceae;Bacillariales;Bacillariaceae;Nitzschia 
                                                                                                       2 
                           Eukaryota;Chromista;Chromobiota;Bacillariophyta;Mediophyceae;Thalassiosirales 
                                                                                                       1 
                                                                                              Unassigned 
                                                                                                       8                                                                                                    