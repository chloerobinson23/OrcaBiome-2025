### OrcaBiome Paper Analysis
### Stacked assembly classes
### March 16, 2025
### Based off Teresita M. Porter, Dec. 23, 2019

# Load required libraries
library(reshape2)  # dcast function
library(ggplot2)   # for visualization
library(dplyr)     # for data wrangling
library(tidyr) # for pivot_longer()
library(gridExtra)  # for arranging plots
library(viridis)

# Read infile (no NTCs)
A <- read.csv(file="OrcaBiome_results_0.8_cut_genera_removed_noNTC_v2.csv", head=TRUE)


# split into flukeprints and BL and run for phyla


# Subset data for ESVs that are only in flukeprints and blanks
FP <- A[A$SampleType == "FP", ]
BL <- A[A$SampleType == "BL", ]


# Pivot to make ESV matrix (pool across versions)
FP.esv <- reshape2::dcast(FP, SampleName ~ Phylum, value.var = "ESVsize", fun.aggregate = sum)
BL.esv <- reshape2::dcast(BL, SampleName ~ Phylum, value.var = "ESVsize", fun.aggregate = sum)


# Convert counts to proportions (normalize per sample)
proportions_table_FP <- FP.esv %>%
  mutate(across(-SampleName, \(x) x / rowSums(FP.esv[, -1]), .names = "{.col}"))

proportions_table_BL <- BL.esv %>%
  mutate(across(-SampleName, \(x) x / rowSums(BL.esv[, -1]), .names = "{.col}"))

# Calculate mean abundance per phylum across all samples
phyla_means_FP <- proportions_table_FP %>%
  summarise(across(-SampleName, mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Phylum", values_to = "Mean_Abundance") %>%
  arrange(desc(Mean_Abundance))

phyla_means_BL <- proportions_table_BL %>%
  summarise(across(-SampleName, mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Phylum", values_to = "Mean_Abundance") %>%
  arrange(desc(Mean_Abundance))


# Identify the top 10 most abundant phyla
top_phyla_FP <- phyla_means_FP$Phylum[1:10]
top_phyla_BL <- phyla_means_BL$Phylum[1:10]

# Convert data to long format
long_data_FP <- proportions_table_FP %>%
  pivot_longer(cols = -SampleName, names_to = "Phylum", values_to = "Proportion") %>%
  mutate(Phylum = ifelse(Phylum %in% top_phyla_FP, Phylum, "Others"))

long_data_BL <- proportions_table_BL %>%
  pivot_longer(cols = -SampleName, names_to = "Phylum", values_to = "Proportion") %>%
  mutate(Phylum = ifelse(Phylum %in% top_phyla_BL, Phylum, "Others"))


# Recalculate proportions so they sum to 1 after grouping "Others"
long_data_FP <- long_data_FP %>%
  group_by(SampleName, Phylum) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop") %>%
  mutate(Proportion = Proportion * 100)  # Convert to percentage

long_data_BL <- long_data_BL %>%
  group_by(SampleName, Phylum) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop") %>%
  mutate(Proportion = Proportion * 100)  # Convert to percentage


# Modify SampleName to show only the part before the second underscore
long_data_FP <- long_data_FP %>%
  mutate(SampleName = sub("(^[^_]+)_.*", "\\1", SampleName))   # Keep only the part before the second underscore

long_data_BL <- long_data_BL %>%
  mutate(SampleName = sub("(^[^_]+)_.*", "\\1", SampleName))   # Keep only the part before the second underscore


# Reorder phylum factor levels, moving "Others" to the end
long_data_FP <- long_data_FP %>%
  mutate(Phylum = factor(Phylum, 
                        levels = c(setdiff(unique(Phylum), "Others"), "Others")))

long_data_BL <- long_data_BL %>%
  mutate(Phylum = factor(Phylum, 
                        levels = c(setdiff(unique(Phylum), "Others"), "Others")))


# Plot the stacked bar chart for flukeprints
p1 <- ggplot(long_data_FP, aes(x = SampleName, y = Proportion, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.7, color = NA) +
  theme_minimal() +
  labs(title = "",
       x = "Flukeprint ID",
       y = "Relative Abundance (%)") +
  scale_fill_viridis(discrete = TRUE)  +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none")

p1


p2 <- ggplot(long_data_BL, aes(x = SampleName, y = Proportion, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.7, color = NA) +
  theme_minimal() +
  labs(title = "",
       x = "Seawater Control ID",
       y = "Relative Abundance (%)") +
  scale_fill_viridis(discrete = TRUE)  +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_blank(), axis.text.y = element_blank())

p2


# Arrange side by side
grid.arrange(p1, p2, ncol = 2)

# Export the arranged plots
ggsave("Fig2_FP_BL_phylum_stacked.png", arrangeGrob(p4, p5, ncol = 2), width = 14, height = 7)


#####
##### Repeat with just flukeprint unique classes
#####


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

ggsave("FigS3_flukeprints_unique_classes.jpeg", plot = p3, dpi = 300, height = 6, width = 8, units = "in")



###### END ######

