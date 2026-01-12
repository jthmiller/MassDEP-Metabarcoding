
## Sample IDs found in the table are missing in the metadata: {'PrimeStrs_7174', 'PrimeStrs_4837'}.


## 2023 data
qiime feature-table filter-features \
    --i-table BB-Tt_MassDiatoms_2023_table.qza \
    --m-metadata-file BB-Tt_MassDiatoms_2023_03012024_hybrid_taxonomy.qza \
    --o-filtered-table filtered_BB-Tt_MassDiatoms_2023_table.qza
######################################


## merge tables for MA 2022 and 2023
qiime feature-table merge \
  --i-tables runs/AlgaeME-rbcLNX031523/qiime_out/AlgaeME-rbcLNX031523_table.qza \
  --i-tables runs/BB-Tt_MassDiatoms_2023/qiime_out/BB-Tt_MassDiatoms_2023_table.qza \
 --o-merged-table results/MA_rbcl_table.qza

## merge asvs for 2022 and 2023
qiime feature-table merge-seqs \
  --i-data runs/AlgaeME-rbcLNX031523/qiime_out/AlgaeME-rbcLNX031523_rep-seqs.qza \
  --i-data runs/BB-Tt_MassDiatoms_2023/qiime_out/BB-Tt_MassDiatoms_2023_rep-seqs.qza \
  --o-merged-data results/MA_rbcl_rep-seqs.qza

qiime feature-table filter-samples \
  --i-table results/MA_rbcl_table.qza \
  --m-metadata-file metadata/MA_2022-2023-metadata.tsv \
  --o-filtered-table results/filtered_MA_rbcl_table.qza

qiime feature-table group \
  --i-table results/filtered_MA_rbcl_table.qza \
  --p-axis sample \
  --m-metadata-file metadata/MA_2022-2023-metadata.tsv \
  --m-metadata-column ID1 \
  --p-mode sum \
  --o-grouped-table results/renamed_MA_rbcl_table.qza


## merge tables for MA 2022,2023,2024
qiime feature-table merge \
  --i-tables runs/101724-MassDEP-Diatoms-2024/qiime_out/101724-MassDEP-Diatoms-2024_table.qza \
  --i-tables results/MA_rbcl_table.qza \
 --o-merged-table results/MassDEP_2022-2024_rbcl_table.qza

## merge asvs for 2022 and 2023
qiime feature-table merge-seqs \
  --i-data runs/101724-MassDEP-Diatoms-2024/qiime_out/101724-MassDEP-Diatoms-2024_rep-seqs.qza \
  --i-data results/MA_rbcl_rep-seqs.qza \
  --o-merged-data results/MassDEP_2022-2024_rbcl_rep-seqs.qza


qiime tools export \
  --input-path results/MassDEP_2022-2024_rbcl_rep-seqs.qza \
  --output-path results/MassDEP_2022-2024_rbcl_rep-seqs


## Do I have this metadata?
qiime feature-table filter-samples \
  --i-table results/MassDEP_2022-2024_rbcl_table.qza \
  --m-metadata-file metadata/MA_2022-2023-renamed-metadata.tsv \
  --o-filtered-table results/filtered_MassDEP_2022-2024_rbcl_table.qza



##############################################################
##############################################################
## merge tables for MA 2022-2025
conda activate qiime2-amplicon-2024.5
cp /home/users/jtm1171/MassDEP/rbcl/results/MassDEP_2022-2024_rbcl_table.qza results/MassDEP_2022-2024_rbcl_table.qza



cut -f1-2 -d',' metadata/rename-2025.csv | tail -n +2 > metadata/rename-2025.tsv
tr ',' '\t' < metadata/rename-2025.csv > metadata/rename-2025.tsv
tr ',' '\t' < metadata/rename-2022-24.csv > metadata/rename-2022-24.tsv

cut -f1-6 -d','  metadata/rename-2022-24.csv | tr ',' '\t'  > metadata/rename-2022-24.tsv


qiime feature-table group \
  --i-table results/MassDEP_2022-2024_rbcl_table.qza \
  --p-axis sample \
  --m-metadata-file metadata/rename-2022-24.tsv \
  --m-metadata-column newname \
  --p-mode sum \
  --o-grouped-table results/MassDEP_2022-2024-renamed_table.qza

qiime feature-table group \
  --i-table runs/MassDiatoms-2025/qiime_out/MassDiatoms-2025_table.qza \
  --p-axis sample \
  --m-metadata-file metadata/rename-2025.tsv \
  --m-metadata-column Newname \
  --p-mode sum \
  --o-grouped-table results/MassDEP_2025-renamed_table.qza

qiime feature-table merge \
  --i-tables results/MassDEP_2025-renamed_table.qza \
  --i-tables results/MassDEP_2022-2024-renamed_table.qza \
 --o-merged-table results/MassDEP_2022-2025_rbcl_table.qza

## merge asvs for 2022 and 2023
qiime feature-table merge-seqs \
  --i-data runs/MassDiatoms-2025/qiime_out/MassDiatoms-2025_rep-seqs.qza \
  --i-data results/MassDEP_2022-2024_rbcl_rep-seqs.qza\
  --o-merged-data results/MassDEP_2022-2025_rep-seqs.qza

qiime tools export \
  --input-path results/MassDEP_2022-2025_rep-seqs.qza \
  --output-path results/MassDEP_2022-2025_rep-seqs

tr ',' '\t' <  metadata/MassDEP_2022-2025-renamed.csv > metadata/MassDEP_2022-2025-renamed.tsv

## Do I have this metadata?
qiime feature-table filter-samples \
  --i-table results/MassDEP_2022-2025_rbcl_table.qza \
  --m-metadata-file metadata/MassDEP_2022-2025-renamed.tsv \
  --o-filtered-table results/filtered_MassDEP_2022-2025_rbcl_table.qza

