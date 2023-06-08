###############################
# RFI extremes
###############################

# Full dataset
sample_metadata<- read_csv("ALL_complete_without_repeatedAnimals_wGroups.csv")

# Filtering animals that are medium efficiency:
ps_genus_RFI <- subset_samples(ps_genus, RFI_groups != "Medium" & RFI_groups != "")
sample_metadata <- subset(sample_metadata, RFI_groups != "Medium" & RFI_groups != "")
ps_genus_RFI_clr <- microbiome::transform(ps_genus_RFI, "clr")
ps_genus_RFI_clr_df <- as.data.frame(otu_table(ps_genus_RFI_clr))

#########################################
## Genus_clr - All days together
#########################################

genus_clr_df <- as.data.frame(otu_table(ps_genus_RFI_clr))
#Preparing the data for PERMANOVA analysis
genus_otu.mat.div <-vegdist(genus_clr_df, method='euclidean')
#Check significance of variables of interest differences
permanova_genus_clr <- adonis2(genus_otu.mat.div ~ sample_metadata$RFI_groups, by="margin", data=genus_clr_df, permutations =  9999, method="euclidean", p.adjust.m='TukeyHSD')
permanova_genus_clr

# Plots for RFI_groups (PCoA and NMDS)
RFI_groups <- betadisper(genus_otu.mat.div, group = sample_metadata$RFI_groups)

# Calculate variation
PCoA1.variation_RFI_groups <- (RFI_groups$eig[1] / (sum(RFI_groups$eig))) * 100
PCoA2.variation_RFI_groups <- (RFI_groups$eig[2] / (sum(RFI_groups$eig))) * 100

# Create data frame for ggplot
RFI_groups_df <- data.frame(PCoA1 = RFI_groups$vectors[, "PCoA1"],
                            PCoA2 = RFI_groups$vectors[, "PCoA2"],
                            RFI_Group = RFI_groups$group)

# Create ggplot object
library(ggplot2)
# Rename the levels of the RFI_Group factor
RFI_groups_df$RFI_Group <- factor(RFI_groups_df$RFI_Group, levels = c("Highest", "Lowest"), labels = c("Most", "Least"))

# Calculate the group averages
group_averages <- aggregate(cbind(PCoA1, PCoA2) ~ RFI_Group, data = RFI_groups_df, FUN = mean)

# Create the ggplot
ggplot(RFI_groups_df, aes(x = PCoA1, y = PCoA2, color = RFI_Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.80, size = 0.4) +
  geom_point(data = group_averages, aes(x = PCoA1, y = PCoA2, color = RFI_Group), size = 6) +
  theme_minimal() +
  labs(title = "",
       subtitle = "PERMANOVA, P < 0.001",
       x = paste0("PCoA1 (", round(RFI_groups$eig[1] / sum(RFI_groups$eig) * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(RFI_groups$eig[2] / sum(RFI_groups$eig) * 100, 2), "%)"),
       color = "RFI Group") +
  scale_color_manual(values = c("Most" = "#033266", "Least" = "#FFD966")) +
  theme(plot.subtitle = element_text(hjust = 0.5))

    