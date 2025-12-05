### Orcabiome analysis
### PICRUSt2 analysis
### December 03, 2025

## Load libraries

library(dplyr)
library(tibble)
library(pheatmap)
library(tidyverse)
library(vegan)
library(readr)
library(viridis)
library(tidyr)

# ========================================
# Filter OTUs unique to flukeprints & sequences
# ========================================

# Load OTU table
otu_table <- read_csv("PiCRUST_OTUs.csv")  # first column = OTU_ID

# Load metadata
metadata <- read_csv("PiCRUST_metadata.csv")  # columns: SampleID, Environment

# Identify flukeprint and seawater samples
flukeprint_ids <- metadata %>%
  filter(Environment == "Flukeprint") %>%
  pull(SampleID)

seawater_ids <- metadata %>%
  filter(Environment == "Seawater") %>%
  pull(SampleID)

# Determine OTUs unique to flukeprints
# Convert sample columns to numeric and replace NAs
otu_table[flukeprint_ids] <- otu_table[flukeprint_ids] %>% mutate(across(everything(), ~ as.numeric(.x))) %>% mutate(across(everything(), ~ replace_na(.x, 0)))
otu_table[seawater_ids] <- otu_table[seawater_ids] %>% mutate(across(everything(), ~ as.numeric(.x))) %>% mutate(across(everything(), ~ replace_na(.x, 0)))

# Sum across flukeprint and seawater samples
otu_table <- otu_table %>%
  rowwise() %>%
  mutate(fluke_sum = sum(c_across(all_of(flukeprint_ids))),
         seawater_sum = sum(c_across(all_of(seawater_ids)))) %>%
  ungroup()

# Keep only OTUs present in flukeprints and absent in seawater
unique_fluke_OTUs <- otu_table %>%
  filter(fluke_sum > 0 & seawater_sum == 0) %>%
  pull(OTU_ID)

# Filter OTU table to unique flukeprint OTUs and only flukeprint columns
otu_fluke_filtered <- otu_table %>%
  filter(OTU_ID %in% unique_fluke_OTUs) %>%
  select(OTU_ID, all_of(flukeprint_ids))

# Save filtered OTU table
write_csv(otu_fluke_filtered, "otu_fluke_unique.csv")

# Filter sequences CSV to match unique OTUs
my_sequences <- read_csv("my_sequences.csv")  # columns: OTU_ID, Sequence

fluke_sequences <- my_sequences %>%
  filter(OTU_ID %in% unique_fluke_OTUs)

write_csv(fluke_sequences, "my_sequences_unique_fluke.csv")


## Convert .csv of ESVs and sequences to FASTA

# Load data
seqs <- read.csv("my_sequences_unique_fluke.csv", stringsAsFactors = FALSE)

# Open a file connection for writing
fasta_file <- file("my_sequences_unique_fluke.fasta", "w")

# Write sequences in FASTA format
for (i in 1:nrow(seqs)) {
  cat(">", seqs$OTU_ID[i], "\n", seqs$Sequence[i], "\n", file = fasta_file, sep = "")
}

# Close connection
close(fasta_file)


# ===============================
# PICRUSt2 Visualization Script
# ===============================

## Load PICRUSt2 output files

# Increase buffer size
Sys.setenv(VROOM_CONNECTION_SIZE = 131072 * 10)  # increase 10×

# EC predictions
bac_ec <- read_tsv("//wsl.localhost/Ubuntu/home/cetusecologist/picrust_Dec_inputs/picrust2_out/bac_EC_predicted.tsv.gz")

# KO predictions
bac_ko <- read_tsv("//wsl.localhost/Ubuntu/home/cetusecologist/picrust_Dec_inputs/picrust2_out/bac_KO_predicted.tsv.gz")

# Pathways
pathways <- read_tsv("//wsl.localhost/Ubuntu/home/cetusecologist/picrust_Dec_inputs/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz")

# NSTI
bac_nsti <- read_tsv("//wsl.localhost/Ubuntu/home/cetusecologist/picrust_Dec_inputs/picrust2_out/bac_marker_predicted_and_nsti.tsv.gz")


# Summarize total abundance per functional feature
ec_totals <- bac_ec %>%
  pivot_longer(-sequence, names_to="EC", values_to="Abundance") %>%
  group_by(EC) %>%
  summarise(Total = sum(Abundance, na.rm=TRUE)) %>%
  arrange(desc(Total))

ko_totals <- bac_ko %>%
  pivot_longer(-sequence, names_to="KO", values_to="Abundance") %>%
  group_by(KO) %>%
  summarise(Total = sum(Abundance, na.rm=TRUE)) %>%
  arrange(desc(Total))

path_totals <- pathways %>%
  pivot_longer(-`pathway`, names_to="Sample", values_to="Abundance") %>%
  group_by(`pathway`) %>%
  summarise(Total = sum(Abundance, na.rm=TRUE)) %>%
  arrange(desc(Total))

# NSTI tells you how closely each OTU matches a reference genome
summary(bac_nsti$metadata_NSTI)
hist(bac_nsti$metadata_NSTI, breaks=30, main="NSTI Distribution", xlab="NSTI")
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.00015  0.02185  0.05809  0.33547  0.12925 26.84915 

# Plot top features
ec_class <- tribble(
  ~EC_prefix, ~Class,
  "1.", "Oxidoreductase",
  "2.", "Transferase",
  "3.", "Hydrolase",
  "4.", "Lyase",
  "5.", "Isomerase",
  "6.", "Ligase"
)


# Select top 20 ECs
top20_ec <- ec_totals %>%
  slice_head(n = 20)

top20_ec


# Add EC class descriptions
top20_ec <- ec_totals %>% slice_head(n=20)


# Define custom colors for EC classes
ec_colors <- c(
  "Oxidoreductase" = "#440154FF",
  "Transferase"     = "#404788FF",
  "Hydrolase"       = "#238A8DFF",
  "Lyase"           = "#29AF7FFF",
  "Isomerase"       = "#95D840FF",
  "Ligase"          = "#FDE725FF"
)

# Prepare top 20 ECs
top20_ec <- ec_totals %>%
  slice_head(n = 20) %>%
  mutate(
    Class = str_extract(EC, "^[0-9]+\\.") %>% 
      str_remove("\\.") %>% 
      as.character() %>% 
      recode(
        "1" = "Oxidoreductase",
        "2" = "Transferase",
        "3" = "Hydrolase",
        "4" = "Lyase",
        "5" = "Isomerase",
        "6" = "Ligase"
      ),
    # Set factor levels to control legend order
    Class = factor(Class, levels = c(
      "Oxidoreductase", 
      "Transferase", 
      "Hydrolase", 
      "Lyase", 
      "Isomerase", 
      "Ligase"
    ))
  )

# Plot bar chart with custom colors
ggplot(top20_ec, aes(x = reorder(EC, Total), y = Total, fill = Class)) +
  geom_col(width = 0.7, color = "black") +          # thicker bars with black borders
  coord_flip() +
  scale_fill_manual(values = ec_colors) +
  labs(
    x = "EC Number",
    y = "Total Predicted Abundance",
    fill = "Enzyme Class",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_blank(),  # remove horizontal grid lines
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "top20_EC_barplot.png",
  plot = last_plot(),
  width = 8,
  height = 6
)


top20_ko <- ko_totals %>% slice_head(n=20)

ggplot(top20_ko, aes(x=reorder(KO, Total), y=Total)) +
  geom_col(fill="darkgreen") +
  coord_flip() +
  labs(x="KO", y="Predicted Abundance", title="Top 20 KEGG Orthologies") +
  theme_minimal()

top20_path <- path_totals %>% slice_head(n=20)

ggplot(top20_path, aes(x=reorder(`pathway`, Total), y=Total)) +
  geom_col(fill="purple") +
  coord_flip() +
  labs(x="Pathway", y="Total Abundance", title="Top 20 Predicted Pathways") +
  theme_minimal()


## Flagging "likely skin functions"

# Use MetaCyc (https://metacyc.org/) to look up pathways and ECs to see if they are associated with skin or not

## Heatmaps of pathway abundance across samples (top 20 only)

# Calculate total abundance for each pathway
top_pathways <- pathways %>%
  mutate(Total = rowSums(across(-pathway), na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  slice_head(n = 20) %>%
  select(-Total)  # remove the total column after selection

# Convert to matrix with pathway names as rownames
path_mat <- top_pathways %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Create a vector matching rownames to annotation
non_skin_pathways <- c(
  "PWY66-389",
  "NONOXIPENT-PWY",
  "PWY-1042"
)

row_annotation <- data.frame(
  SkinAssociated = ifelse(rownames(path_mat) %in% non_skin_pathways, "No", "Yes")
)
rownames(row_annotation) <- rownames(path_mat)


# Define custom colors for Yes/No
ann_colors <- list(
  SkinAssociated = c(
    "No" = turbo(7)[6],  # pick a greenish shade
    "Yes"  = turbo(7)[4]   # pick a reddish shade
  )
)

# Define custom color ramp
my_colors <- colorRampPalette(c("#FDE725FF", "#2A788EFF", "#440154FF"))(50)

dev.off()

# Save as high-resolution PNG
png("heatmap_skin_pathways.png", 
    width = 2000,    # pixels
    height = 2500,   # pixels
    res = 300)       # DPI for publication quality

# Plot heatmap
pheatmap(log1p(path_mat),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         fontsize_row = 8,
         fontsize_col = 8,
         color = my_colors)

dev.off()


## Stats

# Sum abundances across all skin-associated pathways and general pathways per sample
agg <- pathways %>%
  pivot_longer(cols = -pathway, names_to = "Sample", values_to = "Abundance") %>%
  mutate(SkinAssociated = ifelse(pathway %in% non_skin_pathways, "No", "Yes")) %>%
  group_by(Sample, SkinAssociated) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups="drop")

# Pivot to wide format for paired test
wide <- agg %>%
  pivot_wider(names_from = SkinAssociated, values_from = TotalAbundance)

# Calculate differences: Yes - No per sample
diffs <- wide$No - wide$Yes

# Test for normality
shapiro.test(diffs)

#data:  diffs
#W = 0.91569, p-value = 0.1084

# Paired t-test
t.test(wide$No, wide$Yes, paired = TRUE)

#data:  wide$Yes and wide$No
#t = -5.9707, df = 17, p-value = 1.519e-05
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval: -45897.05 -21929.81
#sample estimates:
#mean difference: -33913.43 


#t = -5.9707 → the magnitude of difference relative to variability is large
#p = 1.519e-05 < 0.05 → reject the null hypothesis
#There is a significant difference between skin-associated and general pathways.
#Mean difference = -33913.43 → on average, skin-associated pathways are over 33,000 units higher than general pathways.


########## END ##########

