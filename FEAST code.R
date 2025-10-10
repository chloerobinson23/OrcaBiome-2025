### Orcabiome FEAST analysis
### September 13th, 2025

# Install packages (including interdependencies)
if(!require(devtools)) install.packages("devtools")
library(devtools)

Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
install.packages(Packages)
lapply(Packages, library, character.only = TRUE)

devtools::install_github("cozygene/FEAST")
library(FEAST)


## Data checks!
# Follow GitHub code and set up: https://github.com/cozygene/FEAST
# otus - need to have samples as columns, genera as columns, delete columns with all 0s
# metadata - need to put sinks first, ensure SampleID matches otu table
# IMPORTANT: Two sinks cant share the same id in metadata (copy and paste to duplicate sources for unique ids)

# Load data
metadata <- Load_metadata(metadata_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/metadata_for_FEAST.txt")
otus <- Load_CountMatrix(CountMatrix_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/data_for_FEAST.txt")

FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/",
                      outfile="flukeprint_results")

## Visualize known vs unknown sources

# 1. Load your FEAST output
FEAST_file <- "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/flukeprint_results_source_contributions_matrix.txt"
mixmat <- read.table(FEAST_file, header = TRUE, sep = "\t", row.names = 1)

# 2. Check the first few rows
head(mixmat)

# 3. If the last column is 'unknown', you can calculate the fraction of each sink coming from unknown
unknown_fraction <- mixmat[ , "Unknown"]

# 4. You can also calculate the fraction coming from known sources
known_fraction <- 1 - unknown_fraction

# 5. Combine into a quick summary table
summary_table <- data.frame(
  Sink = rownames(mixmat),
  Fraction_known_sources = known_fraction,
  Fraction_unknown = unknown_fraction
)

# 6. View it
print(summary_table)


# Save as CSV
output_path <- "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/FEAST_summary_table.txt"
write.csv(summary_table, file = sub("\\.txt$", ".csv", output_path), row.names = TRUE)


#### Repeat with just CA and BC whales

## CA

# Load data
metadata_CA <- Load_metadata(metadata_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/metadata_for_FEAST_CA.txt")
otus_CA <- Load_CountMatrix(CountMatrix_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/data_for_FEAST_CA.txt")

FEAST_output <- FEAST(C = otus_CA, metadata = metadata_CA, different_sources_flag = 0, dir_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/",
                      outfile="CA_flukeprint_results")

## Visualize known vs unknown sources for CA samples only

# 1. Load your FEAST output
FEAST_file_CA <- "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/CA_flukeprint_results_source_contributions_matrix.txt"
mixmat_CA <- read.table(FEAST_file_CA, header = TRUE, sep = "\t", row.names = 1)

# 2. Check the first few rows
head(mixmat_CA)

# 3. If the last column is 'unknown', you can calculate the fraction of each sink coming from unknown
unknown_fraction <- mixmat_CA[ , "Unknown"]

# 4. You can also calculate the fraction coming from known sources
known_fraction <- 1 - unknown_fraction

# 5. Combine into a quick summary table
CA_summary_table <- data.frame(
  Sink = rownames(mixmat_CA),
  Fraction_known_sources = known_fraction,
  Fraction_unknown = unknown_fraction
)

# 6. View it
print(CA_summary_table)

# Save as CSV
output_path <- "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/CA_FEAST_summary_table.txt"
write.csv(CA_summary_table, file = sub("\\.txt$", ".csv", output_path), row.names = TRUE)

# Plot it

library(reshape2)
library(ggplot2)
library(viridis)
library(dplyr)

# 1. Reshape FEAST results into long format
mixmat_long <- melt(as.matrix(mixmat_CA), varnames = c("Sink", "Source"), value.name = "Contribution")

# Example: relabel messy source names
mixmat_long$Source <- factor(mixmat_long$Source,
                             levels = unique(mixmat_long$Source), # preserve order
                             labels = gsub("_", " ", unique(mixmat_long$Source))) 


# 2. Separate known sources and Unknown
sources <- unique(mixmat_long$Source)
known_sources <- setdiff(sources, "Unknown")

# Assign inferno colours to known sources
known_cols <- viridis(length(known_sources), option = "inferno")

# Add grey for Unknown
col_map <- c(setNames(known_cols, known_sources), Unknown = "#BEBEBE")

# Publication-quality plot
p1 <- ggplot(mixmat_long, aes(x = Sink, y = Contribution, fill = Source)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75, colour = "black", size = 0.2) + # thin black outline
  scale_fill_manual(values = col_map, name = "Source") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "",
    x = "Flukeprint (Sink)",
    y = "Proportion of Contribution"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.box = "horizontal"
  )

p1

# Save high-res output
ggsave("C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/Orcabiome_analysis/FEAST/CA_contributions_FEAST.png",
       plot = p1, width = 10, height = 7, dpi = 600)


## BC

# Load data
metadata_BC <- Load_metadata(metadata_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/metadata_for_FEAST_BC.txt")
otus_BC <- Load_CountMatrix(CountMatrix_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/data_for_FEAST_BC.txt")

FEAST_output <- FEAST(C = otus_BC, metadata = metadata_BC, different_sources_flag = 0, dir_path = "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/",
                      outfile="BC_flukeprint_results")

## Visualize known vs unknown sources for CA samples only

# 1. Load your FEAST output
FEAST_file_BC <- "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/BC_flukeprint_results_source_contributions_matrix.txt"
mixmat_BC <- read.table(FEAST_file_BC, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 2. Check the first few rows
head(mixmat_BC)

# 3. If the last column is 'unknown', you can calculate the fraction of each sink coming from unknown
unknown_fraction <- mixmat_BC[ , "Unknown"]

# 4. You can also calculate the fraction coming from known sources
known_fraction <- 1 - unknown_fraction

# 5. Combine into a quick summary table
BC_summary_table <- data.frame(
  Sink = rownames(mixmat_BC),
  Fraction_known_sources = known_fraction,
  Fraction_unknown = unknown_fraction
)

# 6. View it
print(BC_summary_table)

# Save as CSV
output_path <- "C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/OrcaBiome_analysis/FEAST/BC_FEAST_summary_table.txt"
write.csv(BC_summary_table, file = sub("\\.txt$", ".csv", output_path), row.names = TRUE)

# Plot it

# 1. Reshape FEAST results into long format
mixmat_long <- melt(as.matrix(mixmat_BC), varnames = c("Sink", "Source"), value.name = "Contribution")

# 2. Clean Source labels
mixmat_long <- mixmat_long %>%
  mutate(Source_clean = as.character(Source)) %>%
  mutate(Source_clean = gsub("_", " ", Source_clean)) %>%     # underscores → spaces
  mutate(Source_clean = gsub(" 1$", " *", Source_clean)) %>%  # numbered suffixes
  mutate(Source_clean = gsub(" 2$", " +", Source_clean)) %>%
  mutate(Source_clean = gsub(" 3$", " ^", Source_clean)) %>%
  mutate(Source_clean = gsub(" 4$", " ~", Source_clean))

# 3. Factor Source_clean so levels are consistent
mixmat_long$Source_clean <- factor(mixmat_long$Source_clean, levels = unique(mixmat_long$Source_clean))

# 4. Build colour map based on Source_clean
sources <- levels(mixmat_long$Source_clean)
known_sources <- setdiff(sources, "Unknown")
known_cols <- viridis(length(known_sources), option = "inferno")
col_map <- c(setNames(known_cols, known_sources), Unknown = "grey50")

# 5. Publication-quality plot using Source_clean
p2 <- ggplot(mixmat_long, aes(x = Sink, y = Contribution, fill = Source_clean)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = col_map, name = "Source") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "",
    x = "Flukeprint (Sink)",
    y = "Proportion of Contribution"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.box = "horizontal"
  )

p2


# stack plots on top of one another before saving

library(patchwork)

# Function to stack plots vertically, keeping individual legends
stack_plots_vertically <- function(plot_list, heights = NULL) {
  # Default equal heights if not specified
  if (is.null(heights)) {
    heights <- rep(1, length(plot_list))
  }
  
  # Start with the first plot
  stacked <- plot_list[[1]]
  
  # Loop through the rest and stack using patchwork with guides = "collect" set FALSE
  for (i in 2:length(plot_list)) {
    stacked <- stacked / plot_list[[i]] + 
      plot_layout(heights = heights, guides = "keep")
  }
  
  return(stacked)
}

# Apply function to p1 and p2
combined_plot <- stack_plots_vertically(list(p1, p2))

# Show combined plot
combined_plot

# Save high-res output
ggsave("C:/Users/ChloeRobinson/OneDrive - Ocean Wise Conservation Association/OWCA/Ocean Wise/Research Projects/eDNA/OrcaBiome/Orcabiome_analysis/FEAST/Combined_CA_BC_FEAST.png",
       plot = combined_plot, width = 10, height = 14, dpi = 600)


####### END #########
