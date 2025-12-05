### OrcaBiome Paper Analysis
### Stacked assembly classes
### December 02, 2025
### Based off Teresita M. Porter, Dec. 23, 2019

# Load required libraries
library(reshape2)  # dcast function
library(ggplot2)   # for visualization
library(dplyr)     # for data wrangling
library(tidyr) # for pivot_longer()
library(gridExtra)  # for arranging plots
library(viridis)

# Read infile (no NTCs)
A <- read.csv(file="OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", head=TRUE)


flukeprints_classes <- c("Acidobacteria_Gp3", "Acidobacteria_Gp7", "BRC1_genera_incertae_sedis",
                         "Candidatus Cloacamonas","Mollicutes","Spirochaetia","Terrimicrobia",
                         "Thermoleophilia","Woesearchaeota Incertae Sedis AR16", "WPS-2_genera_incertae_sedis")  # Replace with actual classes from flukeprints

# Filter long_data to only include the classes found in flukeprints
FP_filtered <- A %>%
  filter(Class %in% flukeprints_classes)


# Pivot to make ESV matrix (pool across versions)
A.1.esv <- reshape2::dcast(FP_filtered, SampleName ~ Class, value.var = "ESVsize", fun.aggregate = sum)

# Convert counts to proportions (normalize per sample)
proportions_table_FP <- A.1.esv %>%
  mutate(across(-SampleName, \(x) x / rowSums(A.1.esv[, -1]), .names = "{.col}"))

# Calculate mean abundance per class across all samples
class_means <- proportions_table_FP %>%
  summarise(across(-SampleName, mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Class", values_to = "Mean_Abundance") %>%
  arrange(desc(Mean_Abundance))

# Identify the top 10 most abundant classes
top_classes_FP <- class_means$Class[1:10]

# Convert data to long format
long_data_FP <- proportions_table_FP %>%
  pivot_longer(cols = -SampleName, names_to = "Class", values_to = "Proportion") %>%
  mutate(Class = ifelse(Class %in% top_classes_FP, Class, "Others"))

# Recalculate proportions so they sum to 1 after grouping "Others"
long_data_FP <- long_data_FP %>%
  group_by(SampleName, Class) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop") %>%
  mutate(Proportion = Proportion * 100)  # Convert to percentage

# Modify SampleName to show only the part before the second underscore
long_data_FP <- long_data_FP %>%
  mutate(SampleName = sub("(^[^_]+)_.*", "\\1", SampleName))  # Keep only the part before the second underscore

# Reorder Bacterial_Class factor levels, moving "Others" to the end
long_data_FP <- long_data_FP %>%
  mutate(Class = factor(Class, 
                        levels = c(setdiff(unique(Class), "Others"), "Others")))


# Plot the stacked bar chart
p3 <- ggplot(long_data_FP, aes(x = SampleName, y = Proportion, fill = Class)) +
  geom_bar(stat = "identity", width = 0.7, color = NA) +
  theme_minimal() +
  labs(title = "",
       x = "Sample ID",
       y = "Relative Abundance (%)") +
  scale_fill_viridis(discrete = TRUE)  +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p3

ggsave("FigS6_flukeprints_unique_classes.jpeg", plot = p3, dpi = 300, height = 6, width = 8, units = "in")



###### END ######

