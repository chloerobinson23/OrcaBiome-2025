### OrcaBiome Paper Analysis
### Individual heatmaps
### March 15, 2025

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
#data is the table produced after metaworks pipeline (only FP data)
A <- read.csv(file="OrcaBiome_results_0.8_cut_genera_removed_noNTC_v2.csv", head=TRUE)


# Read the file containing the flukeprints families
flukeprints_family_file <- "Family_FP.csv"  # Replace with your actual file path
flukeprints_family <- read.csv(flukeprints_family_file)$Family # Assuming the column is named 'Bacterial_Class'


# Filter long_data to only include the genera found in flukeprints
FP_family <- A %>%
  filter(Family %in% flukeprints_family)


# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(FP_family, do.call(rbind, str_split(FP_family$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:30] <- c("SampleID","SampleType", "WhaleID","WhaleSex","Population")

#rename data set
FP <- A.1

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
FP.family<-reshape2::dcast(FP, SampleName ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(FP.family) <- FP.family$SampleName
FP.family$SampleName <- NULL

# remove columns with only zeros
FP.notnull<-FP.family[,colSums(FP.family) !=0]

# remove rows with only zeros & edit rownames
FP.notnull2<-FP.notnull[rowSums(FP.notnull) !=0,]

# calculate 15th percentile for rrarefy function
FP.percentile<-quantile(rowSums(FP.notnull2), prob=0.15)
FP.percentile
# 15% 
# 15.5 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
FP.mat <- rrarefy(FP.notnull2, sample=FP.percentile)

# Convert to df
FP.df<-data.frame(FP.mat)  

# Get total ESVs per sample
FP.df$sums<-rowSums(FP.df)

# Move rownames to first column
FP.df2 <- data.frame(FP.df)
setDT(FP.df2, keep.rownames = TRUE)[]  # Keep row names as a column (rn)


# Get separate whale ID cols from the 'rn' column
setDT(FP.df2)[, c("SampleID", "SampleType", "WhaleID", "WhaleSex", "Population") := tstrsplit(rn, "_", fixed = TRUE)]


# create factors
FP.df2$Population <- factor(FP.df2$Population,
                                  levels = c("MB", "BC"),
                                  labels = c("California", "British Columbia"))
FP.df2$WhaleID <- factor(FP.df2$WhaleID,
                     levels = c("CA122A2","CA10","CA58A1","CA50B","CA51A2A","CA51A3",
                                "T046B4","T037A2","T037","T068C1","T023D3","T073A2",
                                "T073A1","T101B","T100F","T101A","T060D",
                                "T060E"),
                     labels = c("CA122A2","CA10","CA58A1","CA50B","CA51A2A","CA51A3",
                                "T046B4","T037A2","T037","T068C1","T023D3","T073A2",
                                "T073A1","T101B","T100F","T101A","T060D",
                                "T060E"))


# remove unneeded columns (check this one and why)
FP.df3 <- FP.df2[,-c(1,40:42,44)]


# Read in data
#data is the table produced after metaworks pipeline (only FP data)
FP <- read.csv(file="FP_top50_families.csv", head=TRUE)

# Step 4: Melt for ggplot
FP.df4 <- melt(FP.df3, id=c("WhaleID","Population"))

# Step 5: Create factor for better ordering of genera in the plot
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


##########
##### Repeat and now pull all flukeprint data including taxa in both sample types


# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:30] <- c("SampleID","SampleType", "WhaleID","WhaleSex","Population")


# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.1.esv<-reshape2::dcast(A.1, SampleName ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.1.esv) <- A.1.esv$SampleName
A.1.esv$SampleName <- NULL

# remove columns with only zeros
esv.notnull<-A.1.esv[,colSums(A.1.esv) !=0]

# remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 18618.15 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to df
df<-data.frame(rare.mat)  

# grab just FP rows
FP <- df[grepl("FP", rownames(df)),]

# Get total ESVs per sample
FP$sums<-rowSums(FP)

# Move rownames to first column
FP2 <- data.frame(FP)
setDT(FP2, keep.rownames = TRUE)[]  # Keep row names as a column (rn)


# Get separate whale ID cols from the 'rn' column
setDT(FP2)[, c("SampleID", "SampleType", "WhaleID", "WhaleSex", "Population") := tstrsplit(rn, "_", fixed = TRUE)]


# create factors
FP2$Population <- factor(FP2$Population,
                            levels = c("MB", "BC"),
                            labels = c("California", "British Columbia"))
FP2$WhaleID <- factor(FP2$WhaleID,
                         levels = c("CA122A2","CA10","CA58A1","CA50B","CA51A3A","CA51A3",
                                    "T046B4","T037A2","T037","T068C1","T023D3","T073A2",
                                    "T073A1","T101B","T100F","T101A","T060D",
                                    "T060E"),
                         labels = c("CA122A2","CA10","CA58A1","CA50B","CA51A3A","CA51A3",
                                    "T046B4","T037A2","T037","T068C1","T023D3","T073A2",
                                    "T073A1","T101B","T100F","T101A","T060D",
                                    "T060E"))


# remove unneeded columns (check this one and why)
FP3 <- FP2[,-c(1,224,226,228)]


#export FP3 to select top 50 families
write.csv(FP3, "FP3.csv", row.names = TRUE)


# Read in data
#data is the table produced after metaworks pipeline (only FP data)
FP_newdata <- read.csv(file="FP_top50_families_2.csv", head=TRUE)

# Step 4: Melt for ggplot
FP4 <- melt(FP_newdata, id=c("WhaleID","Population"))

# Step 5: Create factor for better ordering of genera in the plot
FP4$variable <- factor(FP4$variable, levels=rev(unique(FP4$variable)))

# Filter data for each Population to only include WhaleIDs present in that facet
FP4_filtered <- FP4 %>%
  group_by(Population) %>%
  filter(WhaleID %in% WhaleID[!is.na(WhaleID)])


FP4_filtered$Population <- fct_relevel(FP4_filtered$Population, "California")


# Compare richness by site
p2.tmp <- ggplot(FP4_filtered) +
  geom_tile(aes(x=WhaleID, y=variable, fill=value)) +
  ggtitle("") +
  labs(x="Whale ID", y="Bacteria Families") +
  theme(legend.title=element_blank()) +
  scale_fill_gradient(na.value="lightgrey", low="white", high="#3B4D93") +
  theme_bw() + 
  facet_grid(cols=vars(Population), scales = "free_x", space = "free_x") +
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

l2 <- get_legend(p2.tmp)

f2 <- ggplot(FP4_filtered) +
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

f2


ggsave("FigS8_flukeprint_heatmap.jpeg", plot = f2, dpi = 300, height = 32, width = 26, units = "cm")
# based on normalized/rarefied data
# Flukeprint bacterial families (also present in SW)
# Faceted by population
# Lab controls excluded



################## END ####################
