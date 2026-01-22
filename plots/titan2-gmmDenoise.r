conda activate qiime2-amplicon-2024.5 && R


# TITAN2 analysis of phyloseq object
# Threshold Indicator Taxa ANalysis for identifying indicator taxa associated with environmental gradients

require(phyloseq)
require(TITAN2)
require(speedyseq)
require(qiime2R)
library(doParallel)
library(gmmDenoise)


set.seed(12345)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================
options <- c(
  'results/filtered_MassDEP_2022-2025_rbcl_table.qza',
  'results/MA_rbcl_2022-2025_hybrid-taxonomy_all.qza',
  'results/MassDEP_2022-2025_rep-seqs.qza',
  'results/MA_rbcl',
  'metadata/qiime_MassDEP_2022-2025-renamed.tsv',
  'results/MassDEP_2022-2025_rooted-tree.qza',
  'refdbs/diat_barcode_v10_263bp-tax.qza',
  'results/default_MA_rbcl_vsearch.qza'
  )

### taxon free 
ps <- qza_to_phyloseq(
  features = options[1], 
  tree = options[6],
  metadata = options[5]
  )


### Count frequency and depth of each ASV across all samples
asv_readcount <- rowSums(otu_table(ps))
asv_pos <- rowSums(otu_table(ps)>0)

## Transform read counts to log10 scale
logrc <- log10(asv_readcount)

# gmmDenoise to identify number of components (k) for GMM
# Perform the cross-validation
cv <- gmmcv(logrc, epsilon = 1e-02)
autoplot(cv)
## Supports k = 3

# Fit a GMM with k = 3
mod <- gmmem(logrc, k = 3)

# quantile.gmmem() returns, by default, the upper one-sided 95% confidence limit 
# of the second uppermost component as the statistically validated abundance
# threshold value
thresh <- quantile(mod) # equivalent to `thresh <- quantile.gmmem(mod)`
# The fitted GMM with the threshold value
autoplot(mod, vline = thresh)


# Remove ASV with threshold 
ps <- prune_taxa(taxa_sums(ps) > 10^thresh, ps)

# Remove samples with zero reads
ps <- prune_samples(sample_sums(ps) > 0, ps)

# Remove ASV with low frequency. TITAN requires taxa to be present in at least 4 samples
ps <- prune_taxa(rowSums(otu_table(ps) > 0) >= 4, ps)


# ============================================================================

# ============================================================================
# 2. PREPARE DATA FOR TITAN2
# ============================================================================

# Extract OTU table (taxa as rows, samples as columns)
otu_table_titan <- t(otu_table(ps))

# Extract environmental variable(s) for analysis
# Using TP/TN as primary environmental gradient

dat_titan <- sample_data(ps)$TP	
## TP <- sample_data(ps)$TP

# Verify data structure
cat("OTU table dimensions:", dim(otu_table_titan), "\n")
cat("Environmental variable length:", length(dat_titan), "\n")
cat("Environmental variable range:", range(dat_titan, na.rm = TRUE), "\n")


cat("Running TITAN2 analysis...\n")

# Run TITAN2 with parallel processing for faster computation
# n_cores <- detectCores() - 1
n_cores <- 12
registerDoParallel(cores = n_cores)

# Set up cluster for TITAN2
cl <- makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# Run TITAN2 (this may take a few minutes)
titan_result <- titan(txa = otu_table_titan, 
                      env = dat_titan, 
                      minSplt = 5,          # Minimum number of taxa occurrences
                      numPerm = 1000,  # Number of permutations for significance
                      boot = TRUE, 
                      nBoot = 500,
                      imax = TRUE,
                      memory = TRUE
                      )

stopCluster(cl)

cat("TITAN2 analysis complete!\n")

save.image('TITAN-ASVs-gmmDenoise_TP.RData')


# ============================================================================
# 4. EXTRACT AND SUMMARIZE RESULTS
# ============================================================================

# Summary statistics
summary(titan_result)

# Extract taxa with significant associations
titan_summary <- titan_result$sppmax




# Filter for significant indicator taxa (p < 0.05)
# significant_taxa <- filter(titan_summary,'reliability' > 0.95)

significant_taxa <- titan_summary[!titan_summary[,'filter'] == 0,]
##significant_taxa <- titan_summary[titan_summary[,'reliability'] > 0.95,]

cat("Number of significant indicator taxa: ", nrow(significant_taxa), "\n")
cat("Significant taxa:\n")
print(significant_taxa)



plot_sumz_density(titan_result, ribbon = FALSE, points = TRUE)

p <- plot_sumz_density(titan_result,
    ribbon = TRUE, points = FALSE, sumz1 = FALSE, change_points = FALSE,
    xlabel = expression(paste("Surface Water Total Nitrogen ("*mu*"g/l)"))
)
ggsave("plots/TITAN2_TN_analysis_plots.pdf", p, width = 14, height = 20, dpi = 300)


plot_sumz_density(titan_result,
ribbon = TRUE, points = FALSE, sumz1 = FALSE, change_points = FALSE,
xlabel = expression(paste("Surface Water Total Phosphorus ("*mu*"g/l)"))
)

plot_taxa_ridges(titan_result,grid = FALSE)


p <- plot_sumz(titan_result)
ggsave("plots/TITAN2_TN-sumz_analysis_plots.pdf", p, width = 14, height = 20, dpi = 300)


plot_sumz(titan_result, filter = TRUE)


plot_taxa_ridges(titan_result, n_ytaxa = 20     )

plot_taxa_ridges(titan_result, n_ytaxa = 30     )



# Split taxa into increasing and decreasing indicators
incr_taxa <- significant_taxa[ significant_taxa[,'maxgrp'] == 1, ]  # Increasing with gradient
decr_taxa <- significant_taxa[ significant_taxa[,'maxgrp'] == 2, ]  # Decreasing with gradient

cat("Increasing indicator taxa (positive association):", nrow(incr_taxa), "\n")
cat("Decreasing indicator taxa (negative association):", nrow(decr_taxa), "\n")
