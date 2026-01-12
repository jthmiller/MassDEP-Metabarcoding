classify_23-24.sh





qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/MassDEP_2022-2024_rbcl_rep-seqs.qza  \
  --i-classifier ../../refdbs/rbcl/diat_barcode_v10_263bp-sklearn-classifier_1.4.2.qza \
  --i-reference-reads ../../refdbs/rbcl/diat_barcode_v10_263bp-seqs.qza \
  --i-reference-taxonomy ../../refdbs/rbcl/diat_barcode_v10_263bp-tax.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification results/MassDEP_2022-2024_rbcl_hybrid-taxonomy_all.qza

### Make output tables
export LD_LIBRARY_PATH='/usr/lib/jvm/java-11-openjdk-11.0.25.0.9-3.el9.x86_64/lib/server:$LD_LIBRARY_PATH'
~/code/qiime2_output_tables.r \
  results/MassDEP_2022-2024_rbcl_table.qza \
  results/MassDEP_2022-2024_rbcl_hybrid-taxonomy_all.qza \
  results/MassDEP_2022-2024_rbcl_rep-seqs.qza \
  MA_rbcl \
  metadata/MA_2022-2023-2024-metadata.tsv


