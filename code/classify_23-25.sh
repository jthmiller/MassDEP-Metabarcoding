classify.sh


### sklearn hybrid taxonomy
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/MA_rbcl_rep-seqs.qza  \
  --i-classifier /home/unhAW/jtmiller/watts/ref-database/rbcl/diat_barcode_v10_263bp-sklearn-classifier.qza \
  --i-reference-reads /home/unhAW/jtmiller/watts/ref-database/rbcl/diat_barcode_v10_263bp-seqs.qza \
  --i-reference-taxonomy /home/unhAW/jtmiller/watts/ref-database/rbcl/diat_barcode_v10_263bp-tax.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification results/MA_rbcl_hybrid-taxonomy.qza

rm -fR results/MA_rbcl_hybrid-taxonomy

qiime tools export \
  --input-path results/MA_rbcl_hybrid-taxonomy.qza \
  --output-path results/MA_rbcl_hybrid-taxonomy

### sklearn hybrid taxonomy
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/MA_rbcl_rep-seqs.qza  \
  --i-classifier /home/unhAW/jtmiller/watts/ref-database/rbcl/diat_barcode_v10_263bp-sklearn-classifier.qza \
  --i-reference-reads /home/unhAW/jtmiller/watts/ref-database/rbcl/diat_barcode_v10_263bp-seqs.qza \
  --i-reference-taxonomy /home/unhAW/jtmiller/watts/ref-database/rbcl/diat_barcode_v10_263bp-tax.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification results/MA_rbcl_hybrid-taxonomy_all.qza

rm -fR results/MA_rbcl_hybrid-taxonomy_all

qiime tools export \
  --input-path results/MA_rbcl_hybrid-taxonomy_all.qza \
  --output-path results/MA_rbcl_hybrid-taxonomy_all

grep '0003e2dc8d6e94c1231839119a79222f'  results/MA_rbcl_hybrid-taxonomy_all/*
grep '0003e2dc8d6e94c1231839119a79222f'  results/MA_rbcl_hybrid-taxonomy/*

## https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0470-z
## indicates that 0.7 confidence is best on mock community

qiime feature-table filter-features \
    --i-table results/renamed_MA_rbcl_table.qza \
    --m-metadata-file results/MA_rbcl_hybrid-taxonomy.qza \
    --o-filtered-table results/filtered_renamed_MA_rbcl_table.qza

qiime taxa barplot \
     --i-table results/filtered_renamed_MA_rbcl_table.qza \
     --i-taxonomy results/MA_rbcl_hybrid-taxonomy.qza \
     --m-metadata-file metadata/MA_2022-2023-renamed-metadata.tsv \
     --o-visualization results/MA_2022-2023-renamed_hybrid_taxa-barplot.qza 

qiime2_output_tables.r results/renamed_MA_rbcl_table.qza \
 results/MA_rbcl_hybrid-taxonomy.qza \
 results/MA_rbcl_rep-seqs.qza \
 results/MA_rbcl_ASV_table.csv \
 results/MA_rbcl_Species_table.csv












