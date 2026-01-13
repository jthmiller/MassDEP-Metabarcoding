classify_2022-2025.sh


classify_2025.sh

######################################
conda activate qiime2-amplicon-2024.5

### sklearn hybrid taxonomy
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/MassDEP_2022-2025_rep-seqs.qza  \
  --i-classifier refdbs/diat_barcode_v10_263bp-sklearn-classifier-2024.5.qza \
  --i-reference-reads refdbs/diat_barcode_v10_263bp-seqs.qza \
  --i-reference-taxonomy refdbs/diat_barcode_v10_263bp-tax.qza \
  --p-threads 12 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification results/MA_rbcl_2022-2025_hybrid-taxonomy_all.qza


### Default params hybrid taxonomy
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query runs/MassDiatoms-2025/qiime_out/MassDiatoms-2025_rep-seqs.qza  \
  --i-classifier refdbs/diat_barcode_v10_263bp-sklearn-classifier-2024.5.qza \
  --i-reference-reads refdbs/diat_barcode_v10_263bp-seqs.qza \
  --i-reference-taxonomy refdbs/diat_barcode_v10_263bp-tax.qza \
  --o-classification results/default_MA_rbcl_2022-2025_hybrid-taxonomy_all.qza


## diat_barcode_v10_263bp_mothur.qza, mothur taxonomy

qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query runs/MassDiatoms-2025/qiime_out/MassDiatoms-2025_rep-seqs.qza  \
  --i-classifier refdbs/diat_barcode_v10_263bp-sklearn-classifier-2024.5.qza \
  --i-reference-reads refdbs/diat_barcode_v10_263bp_qiime.qza \
  --i-reference-taxonomy refdbs/mothur-tax.qza \
  --p-threads 12 \
  --o-classification results/default_MA_rbcl_hybrid-taxonomy_all.qza  

## vsearch default params
qiime feature-classifier classify-consensus-vsearch \
    --i-query results/MassDEP_2022-2025_rep-seqs.qza  \
    --i-reference-reads refdbs/diat_barcode_v10_263bp_qiime.qza \
    --i-reference-taxonomy refdbs/mothur-tax.qza \
    --p-threads 20 \
    --o-classification results/default_MA_rbcl_vsearch \
    --o-search-results results/default_MA_rbcl_vsearch-tophits

## alignment
qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/MassDEP_2022-2025_rep-seqs.qza \
   --o-alignment results/MassDEP_2022-2025_aligned-rep-seqs \
   --o-masked-alignment results/MassDEP_2022-2025_masked-aligned-rep-seqs.qza \
   --o-tree results/MassDEP_2022-2025_unrooted-tree.qza \
   --o-rooted-tree results/MassDEP_2022-2025_rooted-tree.qza \
   --p-n-threads 42 \
   --parallel



### Make output tables

conda activate qiime2-amplicon-2024.5

export LD_LIBRARY_PATH='/usr/lib/jvm/java-17-openjdk-17.0.17.0.10-1.el9.x86_64/lib/server:$LD_LIBRARY_PATH'
code/qiime2_output_tables.r \
  results/filtered_MassDEP_2022-2025_rbcl_table.qza \
  results/MA_rbcl_2022-2025_hybrid-taxonomy_all.qza \
  results/MassDEP_2022-2025_rep-seqs.qza \
  results/MA_rbcl \
  metadata/qiime_MassDEP_2022-2025-renamed.tsv

  