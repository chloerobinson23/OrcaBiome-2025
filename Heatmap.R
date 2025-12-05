### OrcaBiome Paper Analysis
### Individual heatmaps
### December 02, 2025

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(gridExtra) # grid.arrange
library(cowplot) # get_legend
library(dplyr) # group function
library(forcats)  # change order of facets

# Read in data
A <- read.csv(file="OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", head=TRUE)

# Read the file containing the flukeprints families
flukeprints_family_file <- "Family_FP.csv"  # Replace with your actual file path
flukeprints_family <- read.csv(flukeprints_family_file)$Family # Assuming the column is named 'Bacterial_Class'

# Filter long_data to only include the genera found in flukeprints
FP_family <- A %>%
  filter(Family %in% flukeprints_family)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(FP_family, do.call(rbind, str_split(FP_family$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID","SampleType", "WhaleID","WhaleSex","Population","Age")

# Rename data set
FP <- A.1

# Pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
FP.family<-reshape2::dcast(FP, SampleName ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# Move sample to rownames then delete
rownames(FP.family) <- FP.family$SampleName
FP.family$SampleName <- NULL

# Remove columns with only zeros
FP.notnull<-FP.family[,colSums(FP.family) !=0]

# Remove rows with only zeros & edit rownames
FP.notnull2<-FP.notnull[rowSums(FP.notnull) !=0,]

# Convert to df
FP.df<-data.frame(FP.notnull2)  

# Get total ESVs per sample
FP.df$sums<-rowSums(FP.df)

# Move rownames to first column
FP.df2 <- data.frame(FP.df)
setDT(FP.df2, keep.rownames = TRUE)[]  # Keep row names as a column (rn)

# Get separate whale ID cols from the 'rn' column
setDT(FP.df2)[, c("SampleID", "SampleType", "WhaleID", "WhaleSex", "Population","Age") := tstrsplit(rn, "_", fixed = TRUE)]

# Create factors
FP.df2$Population <- factor(FP.df2$Population,
                            levels = c("MB", "BC"),
                            labels = c("California", "British Columbia"))
FP.df2$WhaleID <- factor(FP.df2$WhaleID,
                         levels = c("CA122A2","CA10","CA58A1","CA50B","CA51A3A","CA51A3",
                                    "T046B4","T037A2","T037","T068C1","T023D3","T073A2",
                                    "T073A1","T101B","T100F","T101A","T060D",
                                    "T060E"),
                         labels = c("CA122A2","CA10","CA58A1","CA50B","CA51A3A","CA51A3",
                                    "T046B4","T037A2","T037","T068C1","T023D3","T073A2",
                                    "T073A1","T101B","T100F","T101A","T060D",
                                    "T060E"))


# Remove unneeded columns (check this one and why)
FP.df3 <- FP.df2[,-c(1,40:42,44,46)]


# Melt for ggplot
FP.df4 <- melt(FP.df3, id=c("WhaleID","Population"))

# Create factor for better ordering of genera in the plot
FP.df4$variable <- factor(FP.df4$variable, levels=rev(unique(FP.df4$variable)))

# Filter data for each Population to only include WhaleIDs present in that facet
FP.df4_filtered <- FP.df4 %>%
  group_by(Population) %>%
  filter(WhaleID %in% WhaleID[!is.na(WhaleID)])


# Compare richness by site
p.tmp <- ggplot(FP.df4_filtered) +
  geom_tile(aes(x=WhaleID, y=variable, fill=value)) +
  ggtitle("") +
  labs(x="Whale ID", y="Bacteria Families") +
  theme(legend.title=element_blank()) +
  scale_fill_gradient(na.value="lightgrey", low="white", high="#3B4D93") +
  theme_bw() + 
  facet_grid(cols=vars(Population)) +
  guides(fill=guide_legend(title="ESVs")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1.0,vjust = 0, size=7),
        axis.text.y = element_text(size=7),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size=8),
        legend.text = element_text(size=8))

l <- get_legend(p.tmp)

f <- ggplot(FP.df4_filtered) +
  geom_tile(aes(x = WhaleID, y = variable, fill = value)) +
  ggtitle("") +
  labs(x = "WhaleID", y = "Bacteria Families") +
  theme(legend.title = element_blank()) +
  scale_fill_gradient(na.value = "lightgrey", low = "white", high = "#3B4D93") +
  theme_bw() + 
  facet_grid(cols = vars(Population), scales = "free_x", space = "free_x") +
  guides(fill = guide_legend(title = "ESVs")) +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5, size = 12),  # Increased X axis label font size
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),  
    legend.position = "right", 
    strip.placement = "outside",
    strip.background = element_rect(color = "#3B4D93"),
    strip.text.x = element_text(size = 14), 
    strip.text.y = element_text(size = 14),  
  )

f


ggsave("Fig6_flukeprint_heatmap.jpeg", plot = f, dpi = 300, height = 32, width = 26, units = "cm")
# based on normalized/rarefied data
# Unique Flukeprint bacterial families only
# Faceted by population
# Lab controls excluded



################## END ####################

