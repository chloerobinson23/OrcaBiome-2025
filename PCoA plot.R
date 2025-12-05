### OrcaBiome Paper Analysis
### PCoA plot
### December 02, 2025

# Load required libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)    # for Bray-Curtis and PERMANOVA
library(gridExtra)
library(viridis)

# Read data
A <- read.csv(file="OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", head=TRUE)

# split into flukeprints and BL and run for phyla

# Subset data for ESVs that are only in flukeprints and blanks
FP <- A[A$SampleType == "FP", ]
BL <- A[A$SampleType == "BL", ]


# Pivot to make ESV matrix (pool across versions)
FP.esv <- reshape2::dcast(FP, SampleName ~ Phylum, value.var = "ESVsize", fun.aggregate = sum)
BL.esv <- reshape2::dcast(BL, SampleName ~ Phylum, value.var = "ESVsize", fun.aggregate = sum)

# Convert counts to proportions per sample
FP_rel <- FP.esv %>%
  mutate(across(-SampleName, \(x) x / rowSums(FP.esv[, -1])))

BL_rel <- BL.esv %>%
  mutate(across(-SampleName, \(x) x / rowSums(BL.esv[, -1])))

# Combine FP and BL for overall analysis
combined <- bind_rows(FP_rel %>% mutate(SampleType="FP"),
                      BL_rel %>% mutate(SampleType="BL"))

rownames(combined) <- combined$SampleName
combined_matrix <- as.matrix(combined[, -c(1, ncol(combined))])  # remove SampleName and SampleType for matrix

# Replace NAs with 0
combined_matrix[is.na(combined_matrix)] <- 0

# Bray-Curtis distance
bray_dist <- vegdist(combined_matrix, method="bray")

# PCoA
pcoa_res <- cmdscale(bray_dist, eig=TRUE, k=2)
pcoa_df <- data.frame(
  SampleName = rownames(pcoa_res$points),
  PC1 = pcoa_res$points[,1],
  PC2 = pcoa_res$points[,2],
  SampleType = combined$SampleType
)

# Percent variance explained
var_exp <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)

# PERMANOVA
adonis_res <- adonis2(bray_dist ~ SampleType, data = combined)
print(adonis_res)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = bray_dist ~ SampleType, data = combined)
          #Df SumOfSqs      R2      F Pr(>F)
#Model     1  0.10312 0.03461 1.0039   0.36
#Residual 28  2.87596 0.96539              
#Total    29  2.97907 1.00000 


# Ensure SampleType is a factor with correct labels
pcoa_df$SampleType <- factor(pcoa_df$SampleType,
                             levels = c("BL", "FP"),  # original values in your data
                             labels = c("Seawater control", "Flukeprint"))

# Check levels
levels(pcoa_df$SampleType)
# Should return: "Seawater control" "Flukeprint"

# Plot PCoA
p1 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Flukeprint" = "#FDE725",
      "Seawater control"   = "#2A9D8F"
    )
  ) +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Sample Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "right"
  )

p1

ggsave("Fig3_PCoA.jpeg", plot = p1, dpi = 300, height = 6, width = 8, units = "in")

######### END #########
