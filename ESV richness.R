### OrcaBiome Paper Analysis
### ESV richness
### March 10, 2025

#Load libraries
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
#data is the table produced after metaworks pipeline

a <- read.table("OrcaBiome_results_0.8_cut_genera_removed_noNTC_v2.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
a.1 <- data.frame(a, do.call(rbind, str_split(a$SampleName,"_")), stringsAsFactors = FALSE)
names(a.1)[26:30] <- c("SampleID","SampleType","WhaleID","Sex","Population")

#Sample type = blank or flukeprint, whaleID = ID of whale, sex = sex of whale, population = CA/BC
#edit numbers in [] to correspond with what cols have the above attribute headings

# Create new column for dcast
a.1$SampleTypeWhaleIDSexPopulation <- paste(a.1$SampleType, a.1$WhaleID, a.1$Sex, a.1$Population, sep="_")

# pivot to make esv matrix (pool across versions)
A.esv <- reshape2::dcast(a.1, SampleTypeWhaleIDSexPopulation ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.esv) <- A.esv$SampleTypeWhaleIDSexPopulation
A.esv$SampleTypeWhaleIDSexPopulation <- NULL

#remove columns with only zeros
esv.notnull<-A.esv[,colSums(A.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 18618.15

# set random seed
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
setDT(df2)[, paste0("S", 1:4) := tstrsplit(rn, "_")]
colnames(df2)[colnames(df2)=="S1"] <- "SampleType"
colnames(df2)[colnames(df2)=="S2"] <- "WhaleID"
colnames(df2)[colnames(df2)=="S3"] <- "Sex"
colnames(df2)[colnames(df2)=="S4"] <- "Population"

# create factors
df2$SampleType <- factor(df2$SampleType,
                        levels = c("BL", "FP"),
                        labels = c("Seawater Control", "Flukeprint"))
df2$Sex <- factor(df2$Sex,
                         levels = c("M", "F","U","NA"),
                         labels = c("Male", "Female","Unknown","NA"))
df2$Population<- factor(df2$Population,
                      levels = c("MB", "BC"),
                      labels = c("California", "British Columbia"))


# split by population
california <- df2[df2$Population=="California",]
bc <- df2[df2$Population=="British Columbia",]

# Manually select two specific viridis colors
selected_colors <- c("#2A9D8F", "#FDE725") 


# California plot (Box plot with whiskers)
p <- ggplot(california) +
  geom_boxplot(aes(x=california$SampleType, y=california$sums, fill=california$SampleType)) +  # Box plot
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
  scale_fill_manual(values = selected_colors)  # Optional: Add color using viridis for fill

# British Columbia plot (Box plot with whiskers)
s <- ggplot(bc) +
  geom_boxplot(aes(x=bc$SampleType, y=bc$sums, fill=bc$SampleType)) +  # Box plot
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
  scale_fill_manual(values = selected_colors)  # Optional: Add color using viridis for fill

# Arrange the two plots side by side
g <- grid.arrange(p, s, ncol = 2)

# Print the arranged plots
g

ggsave("Fig1_ESV_richness.jpeg", plot = g, dpi = 300, height = 6, width = 8, units = "in")
# based on normalized data

###### END #######