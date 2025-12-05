### OrcaBiome Paper Analysis
### FEAST and PiCRUST data preparation
### December 02, 2025

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

#data

A <- read.table("OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1 <- data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID", "Type","WhaleID","Sex","Population","Age")

#Sample type = blank or flukeprint, whaleID = ID of whale, sex = sex of whale, population = CA/BC
#edit numbers in [] to correspond with what cols have the above attribute headings

# Create new column for dcast
A.1$SampleTypeWhaleIDSexPopulationAge <- paste(A.1$Type, A.1$WhaleID, A.1$Sex, A.1$Population, A.1$Age, sep="_")

# pivot to make esv matrix (pool across versions)
A.esv <- reshape2::dcast(A.1, SampleTypeWhaleIDSexPopulationAge ~ Genus, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.esv) <- A.esv$SampleTypeWhaleIDSexPopulationAge
A.esv$SampleTypeWhaleIDSexPopulationAge <- NULL

#remove columns with only zeros
esv.notnull<-A.esv[,colSums(A.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 9892.75 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat.abund <- rrarefy(esv.notnull2, sample=esv.percentile)

# Remove columns with all zeros after rarefaction
rare.mat.abund_clean <- rare.mat.abund[, colSums(rare.mat.abund) > 0]

# Remove rows with all zeros (unlikely but safe)
rare.mat.abund_clean <- rare.mat.abund_clean[rowSums(rare.mat.abund_clean) > 0, ]

# Check orientation for FEAST
# FEAST expects rows = features (ESVs/genera), columns = samples
# Currently, your rows are samples, columns are ESVs, so transpose
rare.mat.feast <- t(rare.mat.abund_clean)

# Optional - add rownames as first column for CSV
rare.mat.feast.df <- as.data.frame(rare.mat.feast)
rare.mat.feast.df$ESV <- rownames(rare.mat.feast.df)
rare.mat.feast.df <- rare.mat.feast.df[, c("ESV", setdiff(names(rare.mat.feast.df), "ESV"))]

# Export for FEAST
write.csv(rare.mat.feast.df, "FEAST_input.csv", row.names = FALSE)


###### END ######
