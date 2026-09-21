conda activate qiime2-amplicon-2024.5
export LD_LIBRARY_PATH='/usr/lib/jvm/jre/lib/server:$LD_LIBRARY_PATH'
cd /home/users/jtm1171/spider-webs/12s/runs/102825-PopGen-2025-webs-eDNA/qiime_out
# export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.422.b05-2.el9.x86_64/jre/lib/amd64/server:$LD_LIBRARY_PATH'
R 

# both composiEonal and phylogeneEc
# Unifrac
# Weigthed-Unifrac
# FracEon of the tree specific to either 1 or 2
# FracEon of the diversity specific to 1orto2
# 
# - Jaccard higher than Unifrac
#     - communiEes taxa are disEnct but phylogeneEcally related


### https://genoweb.toulouse.inra.fr/~formation/15_FROGS/6-October2016/FROGS_phyloseq_10102016.pdf

# Reps merged blanks 

#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla
## export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH'
## Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH')

Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/jre/lib/server:$LD_LIBRARY_PATH')


library(tidyverse)
library(qiime2R)
library(taxize)
require(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(rJava)
library(xlsx)
require(speedyseq)
require(openxlsx)
require(rBLAST)
require(RFLPtools)
require(Biostrings)

source('../../../../12s/code/blast-functions.r')


options <- c(
  'A102825-PopGen-2025-webs-eDNA_table.qza',
  'hybrid_taxonomy.qza',
  'A102825-PopGen-2025-webs-eDNA_rep-seqs.qza',
  'A102825-PopGen',
  'spider_metadata.tsv',
  'A102825-PopGen-2025-webs-eDNA_mifish_rooted-tree.qza',
  '/home/users/jtm1171/refdbs/mimammal/12S_MiMammal-U-F/blast_seeds_output/12S_MiMammal-U-F_taxonomy.qza',
  'A102825-PopGen-2025-webs-eDNA_vsearch_taxonomy.qza'
  )


############################################################################################
############################################################################################
############################################################################################
# attach metadata
metadata_out <- read.table(options[5], sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = "@")
metadata_out <- metadata_out[grep('^#',metadata_out[,1], invert = T),]
metadata_out$read_count_warning <- ifelse(metadata_out$non.chimeric < 100, 'warning','ok')
metadata_out$percent_warning <- ifelse(metadata_out$percentage.of.input.non.chimeric < 80, 'warning','ok')


############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################
## plot sequences per sample
physeq <- qza_to_phyloseq(
  features = options[1], 
  taxonomy = options[2],
  tree = options[6],
  metadata = options[5]
  )
# tax_table(physeq)[,'Kingdom']  <- gsub("d__",'k__',tax_table(physeq)[,'Kingdom'])
############################################################################################
############################################################################################
############################################################################################



############################################################################################
############################################################################################
############################################################################################
fish <- read_qza(options[7])$data %>% parse_taxonomy()
fish <- fish$Species
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
## Parse taxonomies

#### Hybrid Taxonomy
tax <- read_qza(options[2])$data
tax$Method <- paste0('hybrid_',tax$Method)
rownames(tax) <- tax$Feature.ID
tax[,'Taxon'] <- gsub("d__",'k__',tax[,'Taxon'])
ptax <- tax %>% parse_taxonomy()



### Vesearch Taxonomy
vtax <- read_qza(options[8])$data
vtax$Method <- 'VSEARCH_MF'
rownames(vtax) <- vtax$Feature.ID
vtax[,'Taxon'] <- gsub("d__",'k__',vtax[,'Taxon'])
vtax <- vtax[which(!vtax$Taxon == 'Unassigned'),]
pvtax <- vtax %>% parse_taxonomy()
colnames(pvtax) <- paste0('VS_',colnames(pvtax))





############################################################################################
############################################################################################
############################################################################################
# filter repseqs by length
seq_min <- 50
seq_max <- 400
threads <- 48

feat_table <- read_qza(options[1])$data
reps <- read_qza(options[3])$data


ASV_length <- names(reps)[which(width(reps) > seq_min & width(reps) < seq_max)]
ASV_length <- unique(c(ASV_length,rownames(ptax),rownames(pvtax),rownames(tax_table(physeq))))
print(paste(length(ASV_length), 'ASVs meet the criteria'))

ASV_length <- ASV_length[which(ASV_length %in% rownames(feat_table))]
print(paste(length(ASV_length), 'ASVs in this data'))

reps <- reps[ASV_length]
ASVs <- names(reps)
############################################################################################
############################################################################################
############################################################################################

pvna <- rowSums(!is.na(pvtax))
pna <- rowSums(!is.na(ptax))

level_count <- sapply(ASVs,function(ASV){
 which.max( c(pna[ASV],pvna[ASV]) )
})

use_pna <- names(which(level_count == 1))
use_pvna <- names(which(level_count == 2))

update_ptax <- pvtax[use_pvna,]
colnames(update_ptax) <- colnames(ptax)
update_ptax$Method <- 'VSEARCH_MF'
update_ptax <- data.frame(rbind(update_ptax, data.frame(ptax[use_pna, ],Method=tax[use_pna,'Method' ])))

sum(ASVs %in% rownames(ptax))
sum(ASVs %in% rownames(pvtax))

############################################################################################
############################################################################################
############################################################################################
## generate blastout
## blast full database
Sys.setenv(BLASTDB = '/home/share/databases/ncbi_nt' )
db <- blast('/home/share/databases/ncbi_nt/nt')

blast_outs <- predict(db, reps, BLAST_args=paste('-max_target_seqs 5 -num_threads',threads), custom_format='qseqid evalue bitscore pident sskingdoms sscinames scomnames',verbose = TRUE, keep_tmp = TRUE)
blast_outs <- fmtBlast(bo_in = blast_outs)
ASVs_blasted <- unique(blast_outs$qseqid) 
euk_asvs <- ASVs_blasted[!ASVs_blasted %in% unique(blast_outs[which(blast_outs$sskingdoms == "Bacteria" ),'qseqid'])]

blast_outs_20 <- predict(db, reps[euk_asvs], BLAST_args=paste('-max_target_seqs 20 -num_threads',threads), custom_format='qseqid evalue bitscore pident sskingdoms sscinames scomnames',verbose = TRUE, keep_tmp = TRUE)
blast_outs_20 <- fmtBlast(bo_in = blast_outs_20)

save.image('after.blast.rsave')

load('after.blast.rsave')

############################################################################################
############################################################################################
############################################################################################


sci <- unique(blast_outs$sscinames)
com <- lapply(sci,function(X, blast_res = blast_outs){
    #unique(blast_res[which(blast_res$sscinames == X),'scomnames'])
    #blast_res[which(blast_res$sscinames == X),'scomnames'][1]
    unique(blast_res[which(blast_res$sscinames == X),'scomnames'])
})
com <- unlist(com)
names(com) <- sci


## add common names to taxonomies
ptax$Common <- com[ptax$Species]
pvtax$VS_Common <- com[pvtax$VS_Species]
update_ptax$Common <- com[update_ptax$Species]

############################################################################################
############################################################################################
############################################################################################
### Get the top evalue for each asv
blast5 <- lapply(ASVs_blasted, bhit_sum, input = blast_outs)
blast5 <- data.frame(do.call(rbind,blast5))
rownames(blast5) <- ASVs_blasted



blast20 <- lapply(euk_asvs, bhit_sum, input = blast_outs_20 )
blast20 <- data.frame(do.call(rbind,blast20))
rownames(blast20) <- euk_asvs

indx <- which(!rownames(blast5) %in% rownames(blast20))
blast20 <- rbind(blast20, blast5[indx,])
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
###### Taxonomy out
tax_out <- data.frame(
    ASVID = ASVs,
    Confidence = NA, 
    Consensus = NA,  
    Method = NA,
    Kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = NA,
    Common = NA,
    row.names = ASVs
)
rownames(tax_out) <- ASVs
tax_out[rownames(update_ptax), colnames(update_ptax) ] <- update_ptax
tax_out[rownames(tax), c('Confidence','Consensus')] <- tax[, c('Confidence','Consensus')] 
ASVs_Taxonomy <- tax_out
################################################

############################################################################################
############################################################################################
############################################################################################
## ASVs out with no blast corrections
ASVs_out <- data.frame(
    ASVID = ASVs,
    Confidence = NA, 
    Consensus = NA,  
    Method = NA,
    Kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = NA,
    Common = NA,
    Total_reads = NA,
    Samples_Pos = NA,
    #Qiime_in_GBIF = NA,
    #BLAST_in_GBIF = NA,
    #BLAST_GBIF_GUESS = NA,
    #BLAST_GBIF_GUESS_in_refdb = NA,
    feat_table[ASVs,],
    seq = as.character(reps[ASVs]),
    check.names = FALSE
    )
rownames(ASVs_out) <- ASVs
ASVs_out[rownames(update_ptax), colnames(update_ptax) ] <- update_ptax
ASVs_out[rownames(tax), c('Confidence','Consensus')] <- tax[, c('Confidence','Consensus')] 

#ASVs_out[, colnames(ptax)] <- ptax[rownames(ASVs_out),]
#ASVs_out[,c('Confidence','Consensus','Method') ] <- tax[rownames(ASVs_out),c('Confidence','Consensus','Method')]
### add these updates to other outputs, and drop 'VS_' columns
#ASVs_out[rownames(update_ptax),colnames(update_ptax)] <- update_ptax
ASVs_out[, 'Total_reads'] <- rowSums(feat_table)[rownames(ASVs_out)]
ASVs_out[, 'Samples_Pos'] <- rowSums(feat_table > 0)[rownames(ASVs_out)]
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
# Table on tree based clustering
tree_otus <- tree_glom(physeq, 0.025)
tree_otus <- prune_taxa(taxa_sums(tree_otus) > 0, tree_otus)

tree_out <- data.frame(
  ASV_ID = rownames(tax_table(tree_otus)),
  tax_table(tree_otus),
  num_samples_pos = rowSums(otu_table(tree_otus)>0),
  total_reads = rowSums(otu_table(tree_otus)),
  otu_table(tree_otus),
  check.names = FALSE
  )
tree_out <- tree_out[,-1]
############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################
# Table on clutsering tips
tip_otus_agnes <- tip_glom(physeq,  h = 0.025, hcfun = cluster::agnes)
tip_otus_agnes <- prune_taxa(taxa_sums(tip_otus_agnes) > 0, tip_otus_agnes)
tip_out_agnes <- data.frame(
  ASV_ID = rownames(tax_table(tip_otus_agnes)),
  tax_table(tip_otus_agnes),
  num_samples_pos = rowSums(otu_table(tip_otus_agnes)>0),
  total_reads = rowSums(otu_table(tip_otus_agnes)),
  otu_table(tip_otus_agnes),
  check.names = FALSE
  )
tip_out_agnes <- tip_out_agnes[,-1]

tip_otus_hclust <- tip_glom(physeq,  h = 0.025, hcfun = hclust)
tip_otus_hclust <- prune_taxa(taxa_sums(tip_otus_hclust) > 0, tip_otus_hclust)
tip_out_hclust <- data.frame(
  ASV_ID = rownames(tax_table(tip_otus_hclust)),
  tax_table(tip_otus_hclust),
  num_samples_pos = rowSums(otu_table(tip_otus_hclust)>0),
  total_reads = rowSums(otu_table(tip_otus_hclust)),
  otu_table(tip_otus_hclust),
  check.names = FALSE
  )
tip_out_hclust <- tip_out_hclust[,-1]

############################################################################################
############################################################################################
############################################################################################

## Blast Output
blast_out <- data.frame(
    ASVID = ASVs,
    Confidence = NA, 
    Consensus = NA,  
    Method = NA,
    Kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = NA,
    Common = NA,
    Total_reads = NA,
    Samples_Pos = NA,
    seq = as.character(reps[ASVs])
    )


blast_out[, colnames(update_ptax)] <- update_ptax[rownames(blast_out),]
blast_out[,c('Confidence','Consensus') ] <- tax[rownames(blast_out),c('Confidence','Consensus')]
blast_out[, 'Total_reads'] <- rowSums(feat_table)[rownames(blast_out)]
blast_out[, 'Samples_Pos'] <- rowSums(feat_table > 0)[rownames(blast_out)]

tax_blast <- blast_out

blast_out <- data.frame(cbind(blast_out, blast20[blast_out$ASVID,],feat_table[blast_out$ASVID,]   ))

##### ID fish
bhits_fish <- rownames(blast_out)[grep('YES',blast_out$is.fish )]
blast_fish <- unique(c(rownames(update_ptax),bhits_fish))
non_fish <- blast_out[!rownames(blast_out) %in% blast_fish, ]
non_fish <- non_fish[!is.na(non_fish$pident),]
getcol <- c("ASVID","Total_reads","Samples_Pos","seq","species_king","species_com","species_sci","evalue","pident",colnames(non_fish)[22:length(non_fish[1,])])
non_fish <- non_fish[,getcol]
non_fish <- non_fish[order(non_fish$species_king,non_fish$species_sci),]

## blast_out <- blast_out[rownames(blast_out) %in% blast_fish, ]

tax_out <- data.frame(cbind(tax_blast, pvtax[tax_blast$ASVID,], blast5[tax_blast$ASVID,]))

## only keep ptax, pvtax, and is.fish results in final blast table
############################################################################################
############################################################################################
############################################################################################

metadata_in = metadata_out
ASVs_in = ASVs_out
nfish = non_fish
blast = blast_out

all <- metadata_in[,'sample.id']
ntc <- grep('NTC', all, value = T)
otc <- grep('Ost', all, value = T)
xb <- grep('XB', all, value = T)
pos <- grep('pos-', all, value = T)
neg <- grep('neg-', all, value = T)
blanks_in <- c(ntc,otc,xb, pos, neg)
samps_in <- all[!all%in% blanks_in]
## ASVs_in <- ASVs_in[rowSums(ASVs_in[,all],na.rm = T)>0,]

## all
otus_out_all <- add_tax(Updated_ASV = ASVs_in, samps = all)
otus_out_all <- otus_out_all[!otus_out_all$Kingdom == 'NA',]
## samples
ASVs_in = ASVs_out
ASVs_in <- ASVs_in[rowSums(ASVs_in[,samps_in],na.rm = T)>0,]
otus_out_samps <- add_tax(Updated_ASV = ASVs_in, samps = samps_in)
otus_out_samps <- otus_out_samps[!otus_out_samps$Kingdom == 'NA',]
## blanks
ASVs_in = ASVs_out
ASVs_in <- ASVs_in[rowSums(ASVs_in[,blanks_in],na.rm = T)>0,]


## otus_out_blanks <- add_tax(Updated_ASV = ASVs_in, samps = blanks_in)

############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################
### OUTPUT XLSX ############################################################################
############################################################################################
out <- list(otus_out_samps, ASVs_out, blast_out)

# out <- list(ASVs_out, otus_out,summary_out, tree_out, tip_out_hclust, tip_out_agnes, blast_out)

tax_cols <- colnames(out[[3]])[5:11]

out <- lapply(out, function(df, cols = tax_cols){  
    orders <- order(apply(df[,cols],1,paste0,collapse=""))

    new_out <- df[orders,]
    #print(head(new_out)[,1:20])
    return(new_out)
})

out <- c(out,list(non_fish),list(metadata_out) )

xlsx_out <- paste0(options[4],'.xlsx')

names(out) <- c("Fish OTUs","ASVs","blast_hits","non_fish","metadata")
openxlsx::write.xlsx(out, file = xlsx_out, rowNames=FALSE)
### OUTPUT XLSX #####
#####################
############################################################################################
############################################################################################
############################################################################################





