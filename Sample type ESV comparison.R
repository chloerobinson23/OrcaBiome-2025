### OrcaBiome Paper Analysis
### FP vs BL per population
### March 11, 2025
### Based off Teresita M. Porter, Jan. 2, 2020

#load libraries
library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(dplyr) # group_by
library(ggrepel) #geom_text_repel

# Read in cat.csv
A <- read.csv(file="OrcaBiome_results_0.8_cut_genera_removed_noNTC_v2.csv", head=TRUE)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1 <- data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:30] <- c("SampleID","SampleType","WhaleID","Sex","Population")


# Create custom field for cast
A.1$OrderGenusGlobalESV <- paste(A.1$Order, A.1$Genus, A.1$GlobalESV, sep=";")

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.1.esv <- reshape2::dcast(A.1, SampleName ~ OrderGenusGlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.1.esv) <- A.1.esv$SampleName
A.1.esv$SampleName <- NULL

# remove columns with only zeros
esv.notnull<-A.1.esv[,colSums(A.1.esv) !=0]

# remove rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]

# calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 18618.15 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# convert to df
df <- as.data.frame(rare.mat)

# Total reads across all samples
total_reads <- rowSums(df)
sum(total_reads)
#519,241

#Get total ESVs (unique ESVs)
total_esvs <- colSums(df > 0)  # Count non-zero occurrences per ESV
length(total_esvs[total_esvs > 0])  # Count unique ESVs detected
#2,217

# convert to presence absence
rare.mat[rare.mat>0] <-1

# convert to df
df <- as.data.frame(rare.mat)

# grab just ANTI and ETOH
FP <- df[grepl("FP", rownames(df)),]
BL <- df[grepl("BL", rownames(df)),]

# Sum ESVs across samples
FP_sums <- colSums(FP)
length(FP_sums[FP_sums > 0])
# 1884

BL_sums <- colSums(BL)
length(BL_sums[BL_sums > 0]) 
# 1556

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
FP_sums[FP_sums>0] <-1
BL_sums[BL_sums>0] <-1

# combine into df
FP_BL <- data.frame(cbind(FP_sums, BL_sums))

# move rownames to first column
setDT(FP_BL, keep.rownames = TRUE)[]
names(FP_BL)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(FP_BL, do.call(rbind, str_split(FP_BL$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumFP=sum(FP_sums),
            sumBL=sum(BL_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "FP vs BL genera.csv", row.names = FALSE)

# Remove any zeroes (only shared taxa now)
t3 <- t3[t3$sumFP > 0 & t3$sumBL > 0, ]


dim(t3)  # Check if dataset is empty
summary(t3)  # Check for NA, Inf, or zero values in sumFP and sumBL
any(is.na(t3$sumFP) | is.na(t3$sumBL))  # Check for NA values
any(t3$sumFP <= 0 | t3$sumBL <= 0)  # Check for non-positive values

# Select top 30 genera by total ESV count
t3_filtered <- t3 %>%
  mutate(total_ESV = sumFP + sumBL) %>%
  arrange(desc(total_ESV)) %>%
  head(30)  # Keep only the top 30 genera


# Compare richness by site
p <- ggplot(t3_filtered, aes(x = sumFP, y = sumBL)) +
  geom_point(aes(color = factor(Order)), position = position_jitter()) +
  labs(x = "Flukeprint ESVs", y = "Seawater Control ESVs", colour = "Order") +
  geom_abline(intercept = 0, slope = 1, linetype = 3) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom") +
  guides(color = guide_legend(ncol = 7))  # Set legend to 7 columns


p

ggsave("Fig3_FP_BL_lineplot.jpeg", plot = p, dpi = 300, height = 16, width = 26, units = "cm")
# based on normalized/rarefied data
# Confidently identified genera only
# top 30 genera only
# NTC controls excluded



###### END #######