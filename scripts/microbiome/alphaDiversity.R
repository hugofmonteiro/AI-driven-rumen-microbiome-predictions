# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

# Full dataset
alphaDiversity <- read_csv("alphaDiversity.csv")

# Descriptive statistics
# Density plots
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
# RFI density per country
RFI_density <- ggplot(data = alphaDiversity, aes(x = RFI, group = Site, fill = Site)) +
  geom_density(adjust = 1.5, alpha = .4) +
  theme_ipsum() +
  ylab("") +
  xlab("Residual feed intake (RFI), kg/d") +
  theme(axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, vjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black')) +
  scale_x_continuous(limits=c(-6,6, by = 2.5))
RFI_density

# MFE density per country
MFE_density <- ggplot(data = alphaDiversity, aes(x = MFE, group = Site, fill = Site)) +
  geom_density(adjust = 1.5, alpha = .4) +
  theme_ipsum() +
  ylab("") +
  xlab("Milk fat efficiency (MFE), g/kg DMI") +
  theme(axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, vjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black')) +
  scale_x_continuous(limits=c(30,130, by = 20))
MFE_density

# MPE density per country
MPE_density <- ggplot(data = alphaDiversity, aes(x = MPE, group = Site, fill = Site)) +
  geom_density(adjust = 1.5, alpha = .4) +
  theme_ipsum() +
  ylab("") +
  xlab("Milk protein efficiency (MPE), g/kg DMI") +
  theme(axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, vjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black')) +
  scale_x_continuous(limits=c(30,75, by = 20))
MPE_density


# Alpha-diversity Plots
library(ggplot2)
library(hrbrthemes) # For theme_ipsum()

#RFI (with all observations)
alphaDiversity_RFI <- subset(alphaDiversity, RFI_image_extremes != "Medium")
alphaDiversity_RFI$RFI_group <- recode(alphaDiversity_RFI$RFI_group,
                                        "Lowest" = "Least",
                                        "Highest" = "Most")

## Chao1 Index
# Perform ANOVA
anova_result <- aov(DMI ~ RFI_groups, data = alphaDiversity_RFI)
anova_summary <- summary(anova_result)
p_value <- round(anova_summary[[1]]["RFI_groups", "Pr(>F)"], 4)

# Function to format p-value
format_p_value <- function(p) {
  if (p > 0.0099) {
    return(sprintf("%.2f", p))
  } else if (p <= 0.0099 & p > 0.00099) {
    return("< 0.01")
  } else {
    return("< 0.001")
  }
}

formatted_p_value <- format_p_value(p_value)

# Create the ggplot object
Chao1 <- alphaDiversity_RFI %>%
  ggplot(aes(x = RFI_groups, y = chao1, fill = RFI_groups)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("#022851", "#FFBF00")) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_ipsum() +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  ggtitle("") +
  xlab("") +
  ylab("Chao 1 Index") +
  labs(subtitle = paste("P", formatted_p_value), hjust = 1.0)
Chao1

## Inverse Simpson Index
# Perform ANOVA
anova_result <- aov(diversity_inverse_simpson ~ RFI_groups, data = alphaDiversity_RFI)
anova_summary <- summary(anova_result)
p_value <- round(anova_summary[[1]]["RFI_groups", "Pr(>F)"], 4)

# Function to format p-value
format_p_value <- function(p) {
  if (p > 0.0099) {
    return(sprintf("%.2f", p))
  } else if (p <= 0.0099 & p > 0.00099) {
    return("< 0.01")
  } else {
    return("< 0.001")
  }
}

formatted_p_value <- format_p_value(p_value)

# Create the ggplot object
InverseSimpson <- alphaDiversity_RFI %>%
  ggplot(aes(x = RFI_groups, y = diversity_inverse_simpson, fill = RFI_groups)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("#022851", "#FFBF00")) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_ipsum() +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  ggtitle("") +
  xlab("") +
  ylab("Inverse Simpson") +
  labs(subtitle = paste("P", formatted_p_value), hjust = 1.0)
InverseSimpson

## Shannon Index
# Perform ANOVA
anova_result <- aov(diversity_shannon ~ RFI_groups, data = alphaDiversity_RFI)
anova_summary <- summary(anova_result)
p_value <- round(anova_summary[[1]]["RFI_groups", "Pr(>F)"], 4)

# Function to format p-value
format_p_value <- function(p) {
  if (p > 0.0099) {
    return(sprintf("%.2f", p))
  } else if (p <= 0.0099 & p > 0.00099) {
    return("< 0.01")
  } else {
    return("< 0.001")
  }
}

formatted_p_value <- format_p_value(p_value)

# Create the ggplot object
Shannon <- alphaDiversity_RFI %>%
  ggplot(aes(x = RFI_groups, y = diversity_shannon, fill = RFI_groups)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("#022851", "#FFBF00")) +
  geom_jitter(color = "black", size = 0.2, alpha = 0.9) +
  theme_ipsum() +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  ggtitle("") +
  xlab("") +
  ylab("Shannon Index") +
  labs(subtitle = paste("P =", formatted_p_value), hjust = 1.0)
Shannon
