classify-2023.sh

### 2023 data
qiime feature-classifier classify-consensus-vsearch \
  --i-query out/rep-seqs.qza \
  --i-reference-reads ref-dbs/diat_barcode_v10_263bp-seq.qza \
  --i-reference-taxonomy ref-dbs/diat_barcode_v10_263bp-tax.qza \
  --p-maxaccepts 10 \
  --p-query-cov 0.8 \
  --p-perc-identity 0.8 \
  --p-threads 12 \
  --p-weak-id 0.50 \
  --o-classification out/vsearch_taxonomy

qiime taxa barplot \
  --i-table out/table.qza \
  --i-taxonomy out/vsearch_taxonomy.qza \
  --o-visualization out/taxa-barplot \
  --m-metadata-file out/dns.qza 

qiime2_output_tables.r qiime_out/BB-Tt_MassDiatoms_2023_table.qza qiime_out/BB-Tt_MassDiatoms_2023_03182024_hybrid_taxonomy.qza qiime_out/BB-Tt_MassDiatoms_2023_rep-seqs.qza qiime_out/BB-Tt_MassDiatoms_2023_03182024_ASV_table.csv qiime_out/BB-Tt_MassDiatoms_2023_03182024_Species_table.csv
qiime2_output_tables.r MA_rbcl_table.qza MA_rbcl_hybrid-taxonomy.qza MA_rbcl_rep-seqs.qza MA_rbcl_ASV_table.csv MA_rbcl_Species_table.csv

qiime taxa barplot \
  --i-table out/table.qza \
  --i-taxonomy BB-Tt_MassDiatoms_2023_03182024_hybrid_taxonomy.qza \
  --o-visualization out/taxa-barplot \
  --m-metadata-file out/dns.qza 



qiime taxa barplot \
     --i-table filtered_BB-Tt_MassDiatoms_2023_table.qza \
     --i-taxonomy BB-Tt_MassDiatoms_2023_03012024_hybrid_taxonomy.qza \
     --o-visualization BB-Tt_MassDiatoms_2023_03012024_hybrid_taxa-barplot \
     --m-metadata-file ../../../metadata/corrected_metadata.tsv

######################################