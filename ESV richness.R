### OrcaBiome Paper Analysis
### ESV richness
### December 02, 2025

# Load libraries
library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(cowplot) # get_legend
library(ggpubr) # normality
library(vegetarian) # calc species equivalents
library(gridExtra) # grid arrange of plots)

#####################################################################
# Look at richness
#####################################################################
# Data
A <- read.table("OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1 <- data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID", "Type","WhaleID","Sex","Population","Age")

# Create new column for dcast
A.1$SampleTypeWhaleIDSexPopulationAge <- paste(A.1$Type, A.1$WhaleID, A.1$Sex, A.1$Population, A.1$Age, sep="_")

# Pivot to make esv matrix (pool across versions)
A.esv <- reshape2::dcast(A.1, SampleTypeWhaleIDSexPopulationAge ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Move sample to rownames then delete
rownames(A.esv) <- A.esv$SampleTypeWhaleIDSexPopulationAge
A.esv$SampleTypeWhaleIDSexPopulationAge <- NULL

# Remove columns with only zeros
esv.notnull<-A.esv[,colSums(A.esv) !=0]

# Remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# Calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 9892.75 

# Set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat.abund <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to presence-absence matrix
rare.mat.abund[rare.mat.abund>0] <-1

# Convert to df
df<-data.frame(rare.mat.abund, stringsAsFactors = FALSE)  

# Get total ESVs per sample
df$sums<-rowSums(df)

# Move rownames to first column
df2<-data.frame(df, stringsAsFactors = FALSE)
setDT(df2, keep.rownames = TRUE)[]

# Get separate substrate and siterep cols
setDT(df2)[, paste0("S", 1:5) := tstrsplit(rn, "_")]
colnames(df2)[colnames(df2)=="S1"] <- "Type"
colnames(df2)[colnames(df2)=="S2"] <- "WhaleID"
colnames(df2)[colnames(df2)=="S3"] <- "Sex"
colnames(df2)[colnames(df2)=="S4"] <- "Population"
colnames(df2)[colnames(df2)=="S5"] <- "Age"

# Create factors
df2$Type <- factor(df2$Type,
                        levels = c("BL", "FP"),
                        labels = c("Seawater Control", "Flukeprint"))
df2$Sex <- factor(df2$Sex,
                         levels = c("M", "F","U","NA"),
                         labels = c("Male", "Female","Unknown","NA"))
df2$Population<- factor(df2$Population,
                      levels = c("MB", "BC"),
                      labels = c("California", "British Columbia"))
df2$Age<- factor(df2$Age,
                        levels = c("SA", "OA", "CA", "NA"),
                        labels = c("Sub-Adult", "Old Adult", "Calf","NA"))


# Split by population
california <- df2[df2$Population=="California",]
bc <- df2[df2$Population=="British Columbia",]

# Manually select two specific viridis colors
selected_colors <- c("#2A9D8F", "#FDE725") 


# California plot (Box plot with whiskers)
p <- ggplot(california) +
  geom_boxplot(aes(x=california$Type, y=california$sums, fill=california$Type)) +  # Box plot
  ggtitle("a)") +
  labs(x="Sample Type", y="Bacteria ESV Richness") +
  theme(legend.title=element_blank()) +
  theme_bw() + 
  facet_grid(cols=vars(Population)) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = selected_colors) +   # Apply green colors to points
  scale_fill_manual(values = selected_colors) +  # Optional: Add color using viridis for fill
  scale_y_continuous(limits = c(0, 600)) 

# British Columbia plot (Box plot with whiskers)
s <- ggplot(bc) +
  geom_boxplot(aes(x=bc$Type, y=bc$sums, fill=bc$Type)) +  # Box plot
  ggtitle("b)") +
  labs(x="Sample Type", y="Bacteria ESV Richness") +
  theme(legend.title=element_blank()) +
  theme_bw() + 
  facet_grid(cols=vars(Population)) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = selected_colors) +   # Apply green colors to points
  scale_fill_manual(values = selected_colors) + # Optional: Add color using viridis for fill
scale_y_continuous(limits = c(0, 600)) 

# Arrange the two plots side by side
g <- grid.arrange(p, s, ncol = 2)
scale_y_continuous(limits = c(0, 600)) 

ggsave("Fig2_ESV_richness.jpeg", plot = g, dpi = 300, height = 6, width = 8, units = "in")
# based on normalized data

###### END #######