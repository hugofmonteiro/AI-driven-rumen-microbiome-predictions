BiocManager::install("ALDEx2", force=TRUE)
remotes::install_github("microbiome/mia")
BiocManager::install("Maaslin2")
BiocManager::install("ANCOMBC")

suppressPackageStartupMessages({ # load packages quietly
  library(mia)
  library(patchwork)
  library(tidySummarizedExperiment)
  library(ALDEx2)
  library(Maaslin2)
  library(MicrobiomeStat)
  library(knitr)
  library(tidyverse)
  library(ANCOMBC)
  library(phyloseq)
  library(microbiome)
  library(tidyverse)
})

###############################
# RFI extremes
###############################

# Filtering animals that are medium efficiency:
ps_genus_RFI <- subset_samples(ps_genus, RFI_groups != "Medium" & RFI_groups != "")
ps_genus_RFI.r <- subset_samples(ps_genus.r, RFI_groups != "Medium" & RFI_groups != "")
ps_genus_RFI_clr <- microbiome::transform(ps_genus_RFI, "clr")
ps_genus_RFI_clr_df <- as.data.frame(otu_table(ps_genus_RFI_clr))

# Suppose you have a phyloseq object called my_phyloseq
# 1. Extract ASV sequence names
asv_names <- rownames(otu_table(ps_genus_RFI))

# 2. Set ASV sequence names as row names of the OTU table
rownames(otu_table(ps_genus_RFI)) <- asv_names

# 3. Create a TreeSummarizedExperiment (TSE) file
ps_genus_RFI_tse <- makeTreeSummarizedExperimentFromPhyloseq(ps_genus_RFI)


###############################
# ALDEx2 - ANOVA-Like Differential Expression tool (Version 2)
###############################

x <- aldex.clr(
  reads = assay(ps_genus_RFI_tse),
  conds = colData(ps_genus_RFI_tse)$RFI_groups, 
  mc.samples = 1000, 
  denom = "all",
  verbose = TRUE
)
x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  hist.plot = TRUE,
  verbose = TRUE)
x_effect <- aldex.effect(x, CI = TRUE, verbose = TRUE)
aldex_out <- data.frame(x_tt, x_effect)
par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05
)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(ps_genus_RFI)), Genus = tax_table(ps_genus_RFI)[, "Genus"])
make_unique_underscore <- function(names) {
  unique_names <- character(length(names))
  for (i in seq_along(names)) {
    current_name <- names[i]
    counter <- 1
    while (current_name %in% unique_names) {
      current_name <- paste0(names[i], "_", counter)
      counter <- counter + 1
    }
    unique_names[i] <- current_name
  }
  return(unique_names)
}

# Merge the data frames
aldex_out_with_names <- merge(aldex_out, genus_info, by.x = "row.names", by.y = "Genus_ID")

# Set the row names of the merged data frame to the unique genus names
unique_genus_names <- make_unique_underscore(aldex_out_with_names$Genus)
rownames(aldex_out_with_names) <- unique_genus_names

# Remove the 'Row.names' and 'Genus' columns
aldex_out_with_names <- aldex_out_with_names[, !(names(aldex_out_with_names) %in% c("Row.names", "Genus"))]

# Filter significant features and display the table
aldex_out_with_names %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the Wilcoxon output rather than tt
  select(we.eBH, wi.eBH, effect, overlap) %>%
  kable()


###############################
# ANCOM-BC (Analysis of compositional microbiomes - with bias correction)
###############################

# Run ANCOMBC function
out = ancombc(phyloseq = ps_genus_RFI.r, formula = "RFI_groups",
              p_adj_method = "holm", zero_cut = 0.99, lib_cut = 0,
              group = "RFI_groups", struc_zero = FALSE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = FALSE,
              alpha = 0.05, global = FALSE)

# Extracting DAA variables
res = out$res
significant_taxa_indices <- which(res$diff_abn$RFI_groupsPositive == TRUE)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(ps_genus_RFI.r)), Genus = tax_table(ps_genus_RFI.r)[, "Genus"])
significant_taxa_info <- genus_info[genus_info$Genus_ID %in% rownames(otu_table(ps_genus_RFI.r)[significant_taxa_indices,]), ]

###############################
# MaAslin2, with bias correction
###############################

asv <- t(assay(ps_genus_RFI_tse))
meta_data <- data.frame(colData(ps_genus_RFI_tse))
fit_data <- Maaslin2(
  asv,
  meta_data,
  output = "RFI_groups",
  transform = "NONE",
  fixed_effects = "RFI_groups",
  # random_effects = c(""), # this has been addressed in the RFI calculation
  reference = "RFI_groups,Positive",  
  normalization = "CLR",
  standardize = FALSE,
  min_prevalence = 0 # make sure prev filtering hasn't been already done
)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(ps_genus_RFI)), Genus = tax_table(ps_genus_RFI)[, "Genus"])

# Filter significant taxa
significant_taxa <- filter(fit_data$results, qval <= 0.05)

# Merge significant_taxa and genus_info to display names
significant_taxa_with_names <- merge(significant_taxa, genus_info, by.x = "feature", by.y = "Genus_ID")

# Remove the 'feature' column, which contains the full sequences
significant_taxa_with_names$feature <- NULL

# Display the table
kable(significant_taxa_with_names)

###############################
# LinDA - linear models for differential abundance analysis
###############################

otu.tab <- as.data.frame(assay(ps_genus_RFI_tse))
meta <- as.data.frame(colData(ps_genus_RFI_tse)) %>% select(RFI_groups)
res_linda <- linda(
  otu.tab,
  meta,
  formula = '~RFI_groups',
  alpha = 0.05,
  prev.filter = 0,
  mean.abund.filter = 0)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(ps_genus_RFI)), Genus = tax_table(ps_genus_RFI)[, "Genus"])

# Convert LinDA output to a data frame and add genus names
linda_output_df <- as.data.frame(res_linda$output)
linda_output_with_names <- merge(linda_output_df, genus_info, by.x = "row.names", by.y = "Genus_ID")

# Make the Genus names unique
unique_genus_names <- make.unique(as.character(linda_output_with_names$Genus))

# Set row names as the unique Genus names
rownames(linda_output_with_names) <- unique_genus_names

# Remove the 'row.names' and 'Genus' columns
linda_output_with_names$Row.names <- NULL
linda_output_with_names$Genus <- NULL

# Filter significant results and display the table using kable
library(knitr)
significant_linda <- linda_output_with_names[linda_output_with_names$RFI_groupsPositive.reject == TRUE, ]
kable(significant_linda)


##########################################################################
# Clustering the DAA methods

summ <- full_join(
  rownames_to_column(aldex_out, "Genus") %>%
    select(Genus, aldex2 = wi.eBH),
  rownames_to_column(out$res$diff_abn, "Genus") %>%
    select(Genus, ancombc = RFI_groupsPositive),
    by = "Genus") %>%
  full_join(
    select(fit_data$results, Genus = feature, maaslin2 = qval), 
    by = "Genus") %>%
  full_join(
    rownames_to_column(as.data.frame(res_linda$output), "Genus") %>%
      select(Genus, LinDA = RFI_groupsPositive.reject), 
    by = "Genus") %>%
  mutate(
    across(c(aldex2, maaslin2), ~ .x <= 0.05),
    # the following line would be necessary without prevalence filtering 
    # as some methods output NA
    across(-Genus, function(x) ifelse(is.na(x), FALSE, x)),
    score = rowSums(across(c(aldex2, ancombc, maaslin2, LinDA)))
  )

# This is how it looks like:
kable(head(summ))

# how many genera were identified by each method
summarise(summ, across(where(is.logical), sum)) %>%
  kable()

# Check the column names in the 'summ' data frame
colnames(summ)

# Check the column names in the 'tax_table(Genus_RFI)' data frame
colnames(genus_info)

# Replace sequences by taxonomy names
summ_with_name <- summ %>%
  left_join(genus_info, by = c("Genus" = "Genus_ID")) %>%
  select(-Genus) %>%
  rename(Genus = Genus.y)

# which genera are identified by more than the following number of methods
RFI_groups_DAA_tests <- filter(summ_with_name, score > 0) # Here you can set to the minimum you want, but better > 0 so you see all significant taxa and how many times they were signficant
RFI_groups_DAA_tests
RFI_groups_DAA_tests_df <- as.data.frame(RFI_groups_DAA_tests)
write.csv(RFI_groups_DAA_tests, "/Users/RFI_groups_DAA_Genus_tests.csv")

# Load the tidyr library
library(tidyr)

# Combine the Genus column with the method columns
RFI_Genus_BestVariables_DAA <- RFI_groups_DAA_tests %>%
  select(Genus, aldex2, ancombc, maaslin2, LinDA) %>%
  pivot_longer(cols = aldex2:LinDA, names_to = "method", values_to = "selected") %>%
  filter(selected) %>%
  mutate(combined_name = paste(Genus, method, sep = "_")) %>%
  select(combined_name)

# Rename the column to 'Genus_Method'
colnames(RFI_Genus_BestVariables_DAA) <- 'Genus_Method'

# Print the new dataframe
print(RFI_Genus_BestVariables_DAA)

# Separate the Genus_Method column into Variable and Method
RFI_Genus_BestVariables_DAA <- RFI_Genus_BestVariables_DAA %>%
  separate(Genus_Method, into = c("Variable", "Method"), sep = "_")

# Print the updated dataframe
print(RFI_Genus_BestVariables_DAA)
write.csv(RFI_Genus_BestVariables_DAA, "/Users/RFI_DAA_Genus_ForNetworkAnalysis.csv") 
# Rerun for each taxonomy level and combine significant taxa for the python network script