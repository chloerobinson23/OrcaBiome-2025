### OrcaBiome Paper Analysis
### Read and ESV summaries
### December 04, 2025
### Based off Teresita M. Porter, Jan. 2, 2020

#load libraries
library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(dplyr) # group_by
library(ggrepel) #geom_text_repel

# Read in data
A <- read.csv(file="OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", head=TRUE)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1 <- data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID", "Type","WhaleID","Sex","Population","Age")

# Create custom field for cast
A.1$GlobalESV <- paste(A.1$GlobalESV, sep=";")

# Pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.1.esv <- reshape2::dcast(A.1, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Move sample to rownames then delete
rownames(A.1.esv) <- A.1.esv$SampleName
A.1.esv$SampleName <- NULL

# Remove columns with only zeros
esv.notnull<-A.1.esv[,colSums(A.1.esv) !=0]

# Remove rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]

# Total reads across all samples (pre-rarefaction)
total_reads <- rowSums(esv.notnull2)
sum(total_reads)
#599124

# Get total ESVs (unique ESVs)
total_esvs <- colSums(df > 0)  # Count non-zero occurrences per ESV
length(total_esvs[total_esvs > 0])  # Count unique ESVs detected
#2,189

# Calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 9892.75  

# Set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to df
df <- as.data.frame(rare.mat)

# Total reads across all samples (post-rarefaction)
total_reads <- rowSums(df)
sum(total_reads)
#274,749

# Get total ESVs (unique ESVs)
total_esvs <- colSums(df > 0)  # Count non-zero occurrences per ESV
length(total_esvs[total_esvs > 0])  # Count unique ESVs detected
#2,189

# Per-sample sequencing depth (total reads per sample)
seq_depth <- rowSums(df)   # sum all reads per sample

# View min/max sequencing depth
range(seq_depth)
summary(seq_depth)


####### END ########