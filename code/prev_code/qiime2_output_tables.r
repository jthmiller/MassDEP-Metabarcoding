#!/home/users/jtm1171/.conda/envs/qiime2-amplicon-2024.5/bin/Rscript --vanilla
## export LD_LIBRARY_PATH='/usr/lib/jvm/java-17-openjdk-17.0.15.0.6-2.el9.x86_64/lib/server:$LD_LIBRARY_PATH'



Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/java-17-openjdk/lib/server:$LD_LIBRARY_PATH')

require(tidyverse)
require(qiime2R)
require(taxize)
require(phyloseq)
require(microbiome)
require(microbiomeutilities)
require(rJava)
require(xlsx)
require(openxlsx)
## Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH')
Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH')

### use: this_script.R <feature_table.qza> <taxonomy.qza> <rep-seqs.qza> <ASV_outname.csv> <Species_outname.csv>

options <- commandArgs(trailingOnly = TRUE)


# options <- c("results/NERRS_18s_5_20_24_filtered-table.qza",
#   "results/NERRS_18s_5_20_24_hybrid_taxonomy.qza",
#   "results/NERRS_18s_5_20_24_rep-seqs.qza",
#   "results/NERRS_18s_5_20_24",
#   "metadata/18s_sample_metadata_5-20.tsv")


 #options <- c(
 #  "results/filtered_HIDAR-COI_table_wControls.qza",
 #  "results/HIDAR-COI_midori_hybrid_taxonomy.qza",
 #  "results/HIDAR-COI_rep-seqs.qza",
 #  "results/HDAR_COI_midori"
 #  )

##options <- c(
##  'qiime_out/HD24-1-18SNX041524_table.qza',
##   'qiime_out/HD24-1-18SNX041524_06042024_hybrid_taxonomy.qza',
##    'qiime_out/HD24-1-18SNX041524_rep-seqs.qza',
##    'qiime_out/HD24-1-18SNX041524_06042024',
##    'qiime_out/HD24-1-18SNX041524_dns_export/metadata.tsv'
##)

#options <- c(
#  'results/MBON_TERNS_JUN17_table_filtered.qza',
#  'results/MBON_TERNS_JUN17_hybrid_taxonomy.qza',
#  'results/MBON_TERNS_JUN17_rep-seqs.qza',
#  'MBON_TERNS_JUN17',
#  'metadata/MBON_TERNS_JUN17_fitler-metadata.tsv')
#
#options <- c(
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_table.qza',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_06212024_hybrid_taxonomy.qza',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_rep-seqs.qza',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_dns_export/metadata.tsv') 
#options <- c(
#  'HD23-RERUN-MFNX051424_table.qza',
#  'HD23-RERUN-MFNX051424_07162024_hybrid_taxonomy.qza',
#  'HD23-RERUN-MFNX051424_rep-seqs.qza',
#  'HD23-RERUN-MFNX051424',
#  'HD23-RERUN-MFNX051424_dns_export/metadata.tsv')
#

#options <- c(
#  'results/filtered-by-features_NERRS_18s_5_20_24.qza',
#  'results/NERRS_18s_CRUX_hybrid_taxonomy.qza',
#  'NERRS_18s_5_20_24_rep-seqs.qza',
#  'NERRS_18s',
#  'metadata/qiime-swmp-sample-metadata.tsv')


# # learn hybrid silva
# ps <- qza_to_phyloseq(
#   tree="results/NERRS_18s_5_20_24_rooted-tree.qza", 
#   features='results/filtered-by-features_NERRS_18s_5_20_24.qza', 
#   taxonomy='results/NERRS_18s_CRUX_hybrid_taxonomy.qza',
#   metadata='metadata/swmp-sample-metadata.tsv'
#   )

#options <- c(
#  'results/filtered-by-features_NERRS_18s_5_20_24.qza',
#  'results/NERRS_18s_CRUX_hybrid_taxonomy.qza',
#  'NERRS_18s_5_20_24_rep-seqs.qza',
#  'NERRS_18s',
#  'metadata/qiime-swmp-sample-metadata.tsv')



## plot sequences per sample
physeq <- qza_to_phyloseq(
  features = options[1], 
  taxonomy = options[2],
  )

reps <- read_qza(options[3])$data
tax <- read_qza(options[2])$data
feats <- read_qza(options[1])$data
rownames(tax) <- tax$Feature.ID

## d__ doesn't work
tax_table(physeq)[,'Kingdom'] <- gsub("d__",'',tax_table(physeq)[,'Kingdom'])
tax_table(physeq)[,'Kingdom'] <- gsub("tax=",'',tax_table(physeq)[,'Kingdom'])


### ASV table
ASVs_out <- data.frame(
  ASV_ID = rownames(tax_table(physeq)),
  tax_table(physeq),
  Confidence = tax[rownames(tax_table(physeq)),'Confidence'],
  ASV_Sequence = as.character(reps[rownames(tax_table(physeq))]),
  num_samples_pos = rowSums(otu_table(physeq)>0),
  otu_table(physeq),
  check.names = FALSE
  )

# Species table
## Collapse taxa (should not have sequences- multiple sequences per taxa)
otus <- tax_glom(physeq, taxrank=rank_names(physeq)[7], NArm=F, bad_empty=c(NA, "", " ", "\t"))

otus_out <- data.frame(
  ASV_ID = rownames(tax_table(otus)),
  tax_table(otus),
  num_samples_pos = rowSums(otu_table(otus)>0),
  otu_table(otus),
  check.names = FALSE
  )
otu_file <- paste0(options[4],"_Species_table.csv")
otus_out <- otus_out[,-1]

# read count
reads_out <- data.frame(cbind(colnames(otu_table(physeq)), colSums(otu_table(physeq))))

# tax table
#tax_out <- data.frame(cbind(rownames(tax_table(physeq)), tax_table(physeq)))
tax_out <- tax[rownames(tax_table(physeq)),] %>% parse_taxonomy()

# attach metadata
metadata_out <- read.table(options[5], sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = "@")


ASV_file <- paste0(options[4],"_ASV_table.csv")
otu_file <- paste0(options[4],"_Species_table.csv")
reads_file <- paste0(options[4],"_reads_table.csv")
tax_file <- paste0(options[4],"_tax_table.csv")

write.csv(ASVs_out, file = ASV_file, row.names=FALSE, quote=FALSE)
write.csv(otus_out, file = otu_file, row.names=FALSE, quote=FALSE)
write.csv(reads_out, file = reads_file, row.names=FALSE, quote=FALSE)
write.csv(tax_out, file = tax_file, row.names=FALSE, quote=FALSE)

xlsx_out <- paste0(options[4],'.xlsx')
out <- list(ASVs_out,otus_out,reads_out,tax_out,metadata_out)
names(out) <- c("ASVs","OTUs","readcount","taxonomy","metadata")

openxlsx::write.xlsx(out, file = xlsx_out, rowNames=FALSE)




# write.xlsx(ASVs_out, file = xlsx_out, sheetName="ASVs", row.names=FALSE)
# write.xlsx(otus_out, file = xlsx_out, sheetName="OTUs", append=TRUE, row.names=FALSE)
# write.xlsx(reads_out, file = xlsx_out, sheetName="readcount", append=TRUE, row.names=FALSE)
# write.xlsx(tax_out, file = xlsx_out, sheetName="taxonomy", append=TRUE, row.names=FALSE)

# openxlsx::write.xlsx(otus_out, file = xlsx_out, sheetName="OTUs", append=TRUE, rowNames=FALSE)
# openxlsx::write.xlsx(reads_out, file = xlsx_out, sheetName="readcount", append=TRUE, rowNames=FALSE)
# openxlsx::write.xlsx(tax_out, file = xlsx_out, sheetName="taxonomy", append=TRUE, rowNames=FALSE)