TITAN2_phyloseq.R



# TITAN2 analysis of phyloseq object
# Threshold Indicator Taxa ANalysis for identifying indicator taxa associated with environmental gradients

require(phyloseq)
require(TITAN2)

require(tidyverse)
require(qiime2R)
require(taxize)
require(phyloseq)
require(microbiome)
require(microbiomeutilities)
require(rJava)
require(xlsx)
require(speedyseq)
require(openxlsx)
require(rBLAST)
require(RFLPtools)
require(Biostrings)
library(doParallel)



library(ggplot2)
library(dplyr)
library(tidyr)
library(microbiome)
library(vegan)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggforce)
library(scales)
library(viridis)
library(patchwork)
library(data.table)
library(foreach)

library(parallel)
library(tidymodels)
library(broom)
library(purrr)
library(stringr)
library(cowplot)
library(ggpmisc)
library(ggsci)
library(ggtext)
library(ggnewscale)
library(ggbeeswarm)

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

##physeq <- qza_to_phyloseq(
##  features = options[1], 
##  taxonomy = options[2],
##  tree = options[6],
##  metadata = options[5]
##  )

### taxon free 
physeq <- qza_to_phyloseq(
  features = options[1], 
  tree = options[6],
  metadata = options[5]
  )

## Agglomerate to OTU level
tree_otus <- tree_glom(physeq, 0.03)
tree_otus_default <- tree_glom(physeq, 0.05)
## tip_glom()

reps <- read_qza(options[3])$data
## cubar::seq_to_codons(as.character(reps[1]))

ps <- tree_otus

### Filter to only samples with environmental data
##sample_data(ps) <- sample_data(ps) %>%
##    filter(!is.na(Nitrate_mgL) & !is.na(Phosphate_mgL) & 
##           !is.na(Silicate_mgL) & !is.na(Chlorophyll_ugL))
##
# Remove samples with zero reads
ps <- prune_samples(sample_sums(ps) > 0, ps)

# Remove ASV with low reads
ps <- prune_taxa(taxa_sums(ps) > 100, ps)

# Remove ASV with low frequency. TITAN requires taxa to be present in at least 4 samples
ps <- prune_taxa(rowSums(otu_table(ps) > 0) >= 4, ps)

# ============================================================================
# 2. PREPARE DATA FOR TITAN2
# ============================================================================

# Extract OTU table (taxa as rows, samples as columns)
otu_table_titan <- t(otu_table(ps))

# Extract environmental variable(s) for analysis
# Using TP/TN as primary environmental gradient

TN <- sample_data(ps)$TN	
TP <- sample_data(ps)$TP

# Verify data structure
cat("OTU table dimensions:", dim(otu_table_titan), "\n")
cat("Environmental variable length:", length(TN), "\n")
cat("Environmental variable range:", range(TN, na.rm = TRUE), "\n")

# ============================================================================
# 3. RUN TITAN2 ANALYSIS
# ============================================================================

cat("Running TITAN2 analysis...\n")

# Run TITAN2 with parallel processing for faster computation
# n_cores <- detectCores() - 1
 n_cores <- 12
registerDoParallel(cores = n_cores)

# Set up cluster for TITAN2
cl <- makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# Run TITAN2 (this may take a few minutes)
titan_result_TN <- titan(txa = otu_table_titan, 
                      env = TN, 
                      minSplt = 5,          # Minimum number of taxa occurrences
                      numPerm = 250,  # Number of permutations for significance
                      boot = TRUE, 
                      nBoot = 500,
                      imax = TRUE,
                      memory = TRUE
                      )

titan_result_TP <- titan(txa = otu_table_titan, 
                      env = TP, 
                      minSplt = 5,          # Minimum number of taxa occurrences
                      numPerm = 250,  # Number of permutations for significance
                      boot = TRUE, 
                      nBoot = 500,
                      imax = TRUE,
                      memory = TRUE
                      )                      

stopCluster(cl)

cat("TITAN2 analysis complete!\n")

# ============================================================================
# 4. EXTRACT AND SUMMARIZE RESULTS
# ============================================================================

# Summary statistics
summary(titan_result)

# Extract taxa with significant associations
titan_summary <- titan_result$sppmax

# Filter for significant indicator taxa (p < 0.05)
significant_taxa <- titan_summary %>%
    filter(pmax < 0.05)

cat("Number of significant indicator taxa: ", nrow(significant_taxa), "\n")
cat("Significant taxa:\n")
print(significant_taxa)

# Split taxa into increasing and decreasing indicators
incr_taxa <- significant_taxa %>% filter(direc == "A")  # Increasing with gradient
decr_taxa <- significant_taxa %>% filter(direc == "B")  # Decreasing with gradient

cat("Increasing indicator taxa (positive association):", nrow(incr_taxa), "\n")
cat("Decreasing indicator taxa (negative association):", nrow(decr_taxa), "\n")

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================

# Plot 1: TITAN2 Results Bubble Plot
p1 <- ggplot(titan_summary, aes(x = cp, y = imax)) +
    geom_point(aes(size = pmax, color = ifelse(pmax < 0.05, "Significant", "Not Significant")),
               alpha = 0.6) +
    scale_size_continuous(range = c(2, 8), name = "p-value", 
                         labels = scales::label_number(accuracy = 0.01)) +
    scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "#CCCCCC")) +
    labs(title = "TITAN2 Indicator Taxon Analysis",
         x = "Change Point (Chlorophyll μg/L)",
         y = "Indicator Value (iMax)",
         color = "Status") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Plot 2: Increasing vs Decreasing Indicators
p2 <- titan_summary %>%
    filter(pmax < 0.05) %>%
    ggplot(aes(x = reorder(sppmax, cp), y = cp, fill = factor(direc))) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("A" = "#1B9E77", "B" = "#D95F02"),
                     labels = c("A" = "Increasing", "B" = "Decreasing")) +
    coord_flip() +
    labs(title = "Significant Indicator Taxon Change Points",
         x = "Taxa",
         y = "Change Point (Chlorophyll μg/L)",
         fill = "Association") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Plot 3: Indicator Value Distribution
p3 <- titan_summary %>%
    ggplot(aes(x = imax, fill = ifelse(pmax < 0.05, "Significant", "Not Significant"))) +
    geom_histogram(bins = 30, alpha = 0.7) +
    scale_fill_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "#CCCCCC")) +
    labs(title = "Distribution of Indicator Values",
         x = "Indicator Value (iMax)",
         y = "Frequency",
         fill = "Status") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Plot 4: Abundance vs Indicator Value
p4 <- titan_summary %>%
    ggplot(aes(x = sum.abund, y = imax)) +
    geom_point(aes(color = ifelse(pmax < 0.05, "Significant", "Not Significant")),
               size = 3, alpha = 0.6) +
    scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "#CCCCCC")) +
    scale_x_log10() +
    labs(title = "Taxon Abundance vs Indicator Value",
         x = "Total Abundance (log scale)",
         y = "Indicator Value (iMax)",
         color = "Status") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Combine plots
combined_plot <- (p1 + p2) / (p3 + p4)

print(combined_plot)

# Save the combined plot
ggsave("results/TITAN2_analysis_plots.pdf", combined_plot, width = 14, height = 10, dpi = 300)
cat("Plots saved to results/TITAN2_analysis_plots.pdf\n")

# ============================================================================
# 6. DETAILED RESULT TABLES
# ============================================================================

# Create comprehensive results table
titan_results_table <- titan_summary %>%
    mutate(
        Taxa = sppmax,
        ChangePoint = cp,
        IndicatorValue = imax,
        PValue = pmax,
        Direction = ifelse(direc == "A", "Increasing", "Decreasing"),
        TotalAbundance = sum.abund,
        Frequency = occurL + occurU,
        Significant = ifelse(pmax < 0.05, "Yes", "No")
    ) %>%
    select(Taxa, ChangePoint, IndicatorValue, PValue, Direction, 
           TotalAbundance, Frequency, Significant) %>%
    arrange(PValue)

# Save results table
write.csv(titan_results_table, "results/TITAN2_significant_taxa.csv", row.names = FALSE)
cat("Results table saved to results/TITAN2_significant_taxa.csv\n")

# Print top 20 indicator taxa
cat("\n===== TOP 20 INDICATOR TAXA =====\n")
print(head(titan_results_table, 20))

# ============================================================================
# 7. INTEGRATED TAXA INFORMATION
# ============================================================================

# Get taxonomy information from phyloseq object
tax_data <- tax_table(ps) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa")

# Merge TITAN2 results with taxonomy
if (nrow(significant_taxa) > 0) {
    titan_with_tax <- titan_results_table %>%
        left_join(tax_data, by = c("Taxa" = "Taxa"))
    
    write.csv(titan_with_tax, "results/TITAN2_significant_taxa_with_taxonomy.csv", row.names = FALSE)
    cat("Taxonomy-annotated results saved to results/TITAN2_significant_taxa_with_taxonomy.csv\n")
    
    cat("\n===== SIGNIFICANT TAXA WITH TAXONOMY =====\n")
    print(titan_with_tax[, 1:10])
}

# ============================================================================
# 8. SAVE SESSION RESULTS
# ============================================================================

# Save TITAN2 result object
saveRDS(titan_result, "results/TITAN2_analysis_object.rds")
cat("TITAN2 result object saved to results/TITAN2_analysis_object.rds\n")

cat("\n===== TITAN2 ANALYSIS COMPLETE =====\n")