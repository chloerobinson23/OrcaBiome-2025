### OrcaBiome Paper Analysis
### Venn diagram
### March 11, 2025

# Load libraries
library(VennDiagram)
library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(dplyr) # group_by
library(ggrepel) #geom_text_repel

# Read the data from the CSV file (created in "Sample Type ESV comparison" script)
data <- read.csv("data_for_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint") %>%
  pull(Genus) %>%
  unique()

# Step 3: Create the Venn diagram
plot1 <-venn.plot <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c("Seawater Control", "Flukeprint"),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom color palette from viridis
  fill = c("lightblue", "black"),  # Set the colors directly  # We need 2 colors, one for each sample type
  cex = 1.5,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

dev.off()

# Display the plot
grid.draw(plot1)

ggsave("FigS4_venn_diagram.jpeg", plot = plot1, dpi = 300, height = 12, width = 16, units = "in")
# NTC and PC removed
# rarefied genus data from "Sample Type ESV Comparison" script (t3)

####
## Repeat with paired venn data
####

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

# convert to presence absence
rare.mat[rare.mat>0] <-1

# convert to df
df <- as.data.frame(rare.mat)

# convert to presence absence
rare.mat[rare.mat>0] <-1

# convert to df
df <- as.data.frame(rare.mat)

#### Pull matched samples ("OB05","OB07") and export

# grab just ANTI and ETOH
OB05 <- df[grepl("OB05", rownames(df)),]
OB07 <- df[grepl("OB07", rownames(df)),]

# Sum ESVs across samples
OB05_sums <- colSums(OB05)
length(OB05_sums[OB05_sums > 0])
# 242

OB07_sums <- colSums(OB07)
length(OB07_sums[OB07_sums > 0]) 
# 166

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB05_sums[OB05_sums>0] <-1
OB07_sums[OB07_sums>0] <-1

# combine into df
OB05_OB07 <- data.frame(cbind(OB05_sums, OB07_sums))

# move rownames to first column
setDT(OB05_OB07, keep.rownames = TRUE)[]
names(OB05_OB07)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB05_OB07, do.call(rbind, str_split(OB05_OB07$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB05=sum(OB05_sums),
            sumOB07=sum(OB07_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB05_OB07.csv", row.names = FALSE)

##### OB08 & OB09

# grab just ANTI and ETOH
OB08 <- df[grepl("OB08", rownames(df)),]
OB09 <- df[grepl("OB09", rownames(df)),]

# Sum ESVs across samples
OB08_sums <- colSums(OB08)
length(OB08_sums[OB08_sums > 0])
# 217

OB09_sums <- colSums(OB09)
length(OB09_sums[OB09_sums > 0]) 
# 173

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB08_sums[OB08_sums>0] <-1
OB09_sums[OB09_sums>0] <-1

# combine into df
OB08_OB09 <- data.frame(cbind(OB08_sums, OB09_sums))

# move rownames to first column
setDT(OB08_OB09, keep.rownames = TRUE)[]
names(OB08_OB09)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB08_OB09, do.call(rbind, str_split(OB08_OB09$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB08=sum(OB08_sums),
            sumOB09=sum(OB09_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OBO8_OB09.csv", row.names = FALSE)



##### OB11 & OB12

# grab just ANTI and ETOH
OB11 <- df[grepl("OB11", rownames(df)),]
OB12 <- df[grepl("OB12", rownames(df)),]

# Sum ESVs across samples
OB11_sums <- colSums(OB11)
length(OB11_sums[OB11_sums > 0])
# 193

OB12_sums <- colSums(OB12)
length(OB12_sums[OB12_sums > 0]) 
# 130

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB11_sums[OB11_sums>0] <-1
OB12_sums[OB12_sums>0] <-1

# combine into df
OB11_OB12 <- data.frame(cbind(OB11_sums, OB12_sums))

# move rownames to first column
setDT(OB11_OB12, keep.rownames = TRUE)[]
names(OB11_OB12)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB11_OB12, do.call(rbind, str_split(OB11_OB12$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB11=sum(OB11_sums),
            sumOB12=sum(OB12_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB11_OB12.csv", row.names = FALSE)

##### OB13 & OB14

# grab just ANTI and ETOH
OB13 <- df[grepl("OB13", rownames(df)),]
OB14 <- df[grepl("OB14", rownames(df)),]

# Sum ESVs across samples
OB13_sums <- colSums(OB13)
length(OB13_sums[OB13_sums > 0])
# 318

OB14_sums <- colSums(OB14)
length(OB14_sums[OB14_sums > 0]) 
# 245

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB13_sums[OB13_sums>0] <-1
OB14_sums[OB14_sums>0] <-1

# combine into df
OB13_OB14 <- data.frame(cbind(OB13_sums, OB14_sums))

# move rownames to first column
setDT(OB13_OB14, keep.rownames = TRUE)[]
names(OB13_OB14)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB13_OB14, do.call(rbind, str_split(OB13_OB14$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB13=sum(OB13_sums),
            sumOB14=sum(OB14_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB13_OB14.csv", row.names = FALSE)

##### OB18 & OB19

# grab just ANTI and ETOH
OB18 <- df[grepl("OB18", rownames(df)),]
OB19 <- df[grepl("OB19", rownames(df)),]

# Sum ESVs across samples
OB18_sums <- colSums(OB18)
length(OB18_sums[OB18_sums > 0])
# 186

OB19_sums <- colSums(OB19)
length(OB19_sums[OB19_sums > 0]) 
# 129

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB18_sums[OB18_sums>0] <-1
OB19_sums[OB19_sums>0] <-1

# combine into df
OB18_OB19 <- data.frame(cbind(OB18_sums, OB19_sums))

# move rownames to first column
setDT(OB18_OB19, keep.rownames = TRUE)[]
names(OB18_OB19)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB18_OB19, do.call(rbind, str_split(OB18_OB19$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB18=sum(OB18_sums),
            sumOB19=sum(OB19_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB18_OB19.csv", row.names = FALSE)


##### OB21 & OB22

# grab just ANTI and ETOH
OB21 <- df[grepl("OB21", rownames(df)),]
OB22 <- df[grepl("OB22", rownames(df)),]

# Sum ESVs across samples
OB21_sums <- colSums(OB21)
length(OB21_sums[OB21_sums > 0])
# 293

OB22_sums <- colSums(OB22)
length(OB22_sums[OB22_sums > 0]) 
# 236

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB21_sums[OB21_sums>0] <-1
OB22_sums[OB22_sums>0] <-1

# combine into df
OB21_OB22 <- data.frame(cbind(OB21_sums, OB22_sums))

# move rownames to first column
setDT(OB21_OB22, keep.rownames = TRUE)[]
names(OB21_OB22)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB21_OB22, do.call(rbind, str_split(OB21_OB22$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB21=sum(OB21_sums),
            sumOB22=sum(OB22_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB21_OB22.csv", row.names = FALSE)


##### Start of BC samples - OB23, OB25, OB24

# grab just ANTI and ETOH
OB23 <- df[grepl("OB23", rownames(df)),]
OB25 <- df[grepl("OB25", rownames(df)),]
OB24 <- df[grepl("OB24", rownames(df)),]

# Sum ESVs across samples
OB23_sums <- colSums(OB23)
length(OB23_sums[OB23_sums > 0])
# 711

OB25_sums <- colSums(OB25)
length(OB25_sums[OB25_sums > 0]) 
# 544

OB24_sums <- colSums(OB24)
length(OB24_sums[OB24_sums > 0]) 
# 525

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB23_sums[OB23_sums>0] <-1
OB25_sums[OB25_sums>0] <-1
OB24_sums[OB24_sums>0] <-1

# combine into df
OB23_OB25_OB24 <- data.frame(cbind(OB23_sums, OB25_sums, OB24_sums))

# move rownames to first column
setDT(OB23_OB25_OB24, keep.rownames = TRUE)[]
names(OB23_OB25_OB24)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB23_OB25_OB24, do.call(rbind, str_split(OB23_OB25_OB24$Taxon_GlobalESV,";")))
names(t)[5:7] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB23=sum(OB23_sums),
            sumOB25=sum(OB25_sums),
            sumOB24=sum(OB24_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[5:6] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB23_OB25_OB24.csv", row.names = FALSE)


##### OB26 & OB27

# grab just ANTI and ETOH
OB26 <- df[grepl("OB26", rownames(df)),]
OB27 <- df[grepl("OB27", rownames(df)),]

# Sum ESVs across samples
OB26_sums <- colSums(OB26)
length(OB26_sums[OB26_sums > 0])
# 779

OB27_sums <- colSums(OB27)
length(OB27_sums[OB27_sums > 0]) 
# 702

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB26_sums[OB26_sums>0] <-1
OB27_sums[OB27_sums>0] <-1

# combine into df
OB26_OB27 <- data.frame(cbind(OB26_sums, OB27_sums))

# move rownames to first column
setDT(OB26_OB27, keep.rownames = TRUE)[]
names(OB26_OB27)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB26_OB27, do.call(rbind, str_split(OB26_OB27$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB26=sum(OB26_sums),
            sumOB27=sum(OB27_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB26_OB27.csv", row.names = FALSE)


##### OB29 & OB30

# grab just ANTI and ETOH
OB29 <- df[grepl("OB29", rownames(df)),]
OB30 <- df[grepl("OB30", rownames(df)),]

# Sum ESVs across samples
OB29_sums <- colSums(OB29)
length(OB29_sums[OB29_sums > 0])
# 720

OB30_sums <- colSums(OB30)
length(OB30_sums[OB30_sums > 0]) 
# 401

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB29_sums[OB29_sums>0] <-1
OB30_sums[OB30_sums>0] <-1

# combine into df
OB29_OB30 <- data.frame(cbind(OB29_sums, OB30_sums))

# move rownames to first column
setDT(OB29_OB30, keep.rownames = TRUE)[]
names(OB29_OB30)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB29_OB30, do.call(rbind, str_split(OB29_OB30$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB29=sum(OB29_sums),
            sumOB30=sum(OB30_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB29_OB30.csv", row.names = FALSE)


##### OB31, OB32, OB33, OB34

# grab just ANTI and ETOH
OB31 <- df[grepl("OB31", rownames(df)),]
OB32 <- df[grepl("OB32", rownames(df)),]
OB33 <- df[grepl("OB33", rownames(df)),]
OB34 <- df[grepl("OB34", rownames(df)),]

# Sum ESVs across samples
OB31_sums <- colSums(OB31)
length(OB31_sums[OB31_sums > 0])
# 476

OB32_sums <- colSums(OB32)
length(OB32_sums[OB32_sums > 0]) 
# 416

OB33_sums <- colSums(OB33)
length(OB33_sums[OB33_sums > 0]) 
# 414

OB34_sums <- colSums(OB34)
length(OB34_sums[OB34_sums > 0]) 
# 60

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB31_sums[OB31_sums>0] <-1
OB32_sums[OB32_sums>0] <-1
OB33_sums[OB33_sums>0] <-1
OB34_sums[OB34_sums>0] <-1

# combine into df
OB31_OB32_OB33_OB34 <- data.frame(cbind(OB31_sums, OB32_sums, OB33_sums, OB34_sums))

# move rownames to first column
setDT(OB31_OB32_OB33_OB34, keep.rownames = TRUE)[]
names(OB31_OB32_OB33_OB34)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB31_OB32_OB33_OB34 , do.call(rbind, str_split(OB31_OB32_OB33_OB34 $Taxon_GlobalESV,";")))
names(t)[6:8] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB31=sum(OB31_sums),
            sumOB32=sum(OB32_sums),
            sumOB33=sum(OB33_sums),
            sumOB34=sum(OB34_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[6:7] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB31_OB32_OB33_OB34.csv", row.names = FALSE)


##### OB35, OB36, OB37, OB38, OB39

# grab just ANTI and ETOH
OB35 <- df[grepl("OB35", rownames(df)),]
OB37 <- df[grepl("OB37", rownames(df)),]
OB38 <- df[grepl("OB38", rownames(df)),]
OB39 <- df[grepl("OB39", rownames(df)),]

# Sum ESVs across samples
OB35_sums <- colSums(OB35)
length(OB35_sums[OB35_sums > 0])
# 875

OB37_sums <- colSums(OB37)
length(OB37_sums[OB37_sums > 0]) 
# 605

OB38_sums <- colSums(OB38)
length(OB38_sums[OB38_sums > 0]) 
# 667

OB39_sums <- colSums(OB39)
length(OB39_sums[OB39_sums > 0]) 
# 632

# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB35_sums[OB35_sums>0] <-1
OB37_sums[OB37_sums>0] <-1
OB38_sums[OB38_sums>0] <-1
OB39_sums[OB39_sums>0] <-1

# combine into df
OB35_OB37_OB38_OB39 <- data.frame(cbind(OB35_sums, OB37_sums, OB38_sums, OB39_sums))

# move rownames to first column
setDT(OB35_OB37_OB38_OB39, keep.rownames = TRUE)[]
names(OB35_OB37_OB38_OB39)[1] <- "Taxon_GlobalESV"

# remove last two fields of Taxon
t <- data.frame(OB35_OB37_OB38_OB39 , do.call(rbind, str_split(OB35_OB37_OB38_OB39$Taxon_GlobalESV,";")))
names(t)[6:8] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB35=sum(OB35_sums),
            sumOB37=sum(OB37_sums),
            sumOB38=sum(OB38_sums),
            sumOB39=sum(OB39_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[6:7] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB35_OB37_OB38_OB39.csv", row.names = FALSE)


##### OB40, OB41, OB42

# grab just ANTI and ETOH
OB40 <- df[grepl("OB40", rownames(df)),]
OB41 <- df[grepl("OB41", rownames(df)),]
OB42 <- df[grepl("OB42", rownames(df)),]


# Sum ESVs across samples
OB40_sums <- colSums(OB40)
length(OB40_sums[OB40_sums > 0])
# 793

OB41_sums <- colSums(OB41)
length(OB41_sums[OB41_sums > 0]) 
# 708

OB42_sums <- colSums(OB42)
length(OB42_sums[OB42_sums > 0]) 
# 569


# Convert to presence absence (i.e. pool across samples above, then mark when detected with a 1)
OB40_sums[OB40_sums>0] <-1
OB41_sums[OB41_sums>0] <-1
OB42_sums[OB42_sums>0] <-1


# combine into df
OB40_OB41_OB42 <- data.frame(cbind(OB40_sums, OB41_sums, OB42_sums))

# move rownames to first column
setDT(OB40_OB41_OB42, keep.rownames = TRUE)[]
names(OB40_OB41_OB42)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(OB40_OB41_OB42, do.call(rbind, str_split(OB40_OB41_OB42$Taxon_GlobalESV,";")))
names(t)[5:7] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumOB40=sum(OB40_sums),
            sumOB41=sum(OB41_sums),
            sumOB42=sum(OB42_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[5:6] <- c("Order","Genus")


# Save to .csv (before removing zeroes)
write.csv(t3, "OB40_OB41_OB42.csv", row.names = FALSE)



######## Create paired venn diagrams

library(VennDiagram)
library(dplyr)

### OB05 & OB07
# Step 1: Read the data from the CSV file
data <- read.csv("OB05_OB07_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB07") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB05") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB05-0B07_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot1 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot1)  # Display the Venn diagram



# Close the graphics device and save the file
dev.off()  

# NTC and PC removed
# rarefied genus data from "Sample Type ESV Comparison" script (t3)


### OB08 & OB09
# Step 1: Read the data from the CSV file
data <- read.csv("OB08_OB09_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB09") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB08") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB08-0B09_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot2 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot2)  # Display the Venn diagram


# Close the graphics device and save the file
dev.off()  


### OB11 & OB12

# Step 1: Read the data from the CSV file
data <- read.csv("OB11_OB12_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB12") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB11") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB11-0B12_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot3 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot3)  # Display the Venn diagram


# Close the graphics device and save the file
dev.off()  


### OB13 & OB14

# Step 1: Read the data from the CSV file
data <- read.csv("OB13_OB14_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB14") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB13") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB13-0B14_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot4 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot4)  # Display the Venn diagram


# Close the graphics device and save the file
dev.off()


### OB18 & OB19

# Step 1: Read the data from the CSV file
data <- read.csv("OB18_OB19_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB19") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB18") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB18-0B19_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot5 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot5)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()


### OB21 & OB22

# Step 1: Read the data from the CSV file
data <- read.csv("OB21_OB22_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB22") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB21") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB21-0B22_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot6 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot6)  # Display the Venn diagram


# Close the graphics device and save the file
dev.off()


### OB23 & OB24

# Step 1: Read the data from the CSV file
data <- read.csv("OB23_OB24_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB24") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB23") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB23-0B24_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot7 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot7)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()

### OB25 & OB24

# Step 1: Read the data from the CSV file
data <- read.csv("OB25_OB24_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB24") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB25") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB25-0B24_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot8 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot8)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()

### OB26 & OB27

# Step 1: Read the data from the CSV file
data <- read.csv("OB26_OB27_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB27") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB26") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB26-0B27_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot9 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot9)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()


### OB29 & OB30

# Step 1: Read the data from the CSV file
data <- read.csv("OB29_OB30_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB30") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB29") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB29-0B30_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot10 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot10)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()



### OB31 & OB34

# Step 1: Read the data from the CSV file
data <- read.csv("OB31_OB34_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB34") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB31") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB31-OB34_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot11 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot11)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()


### OB32 & OB34

# Step 1: Read the data from the CSV file
data <- read.csv("OB32_OB34_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB34") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB32") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB32-OB34_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot12 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot12)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()



### OB33 & OB34

# Step 1: Read the data from the CSV file
data <- read.csv("OB33_OB34_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB34") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB33") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB33-OB34_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot13 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot13)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()




### OB35 & OB39

# Step 1: Read the data from the CSV file
data <- read.csv("OB35_OB39_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB39") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB35") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB35-OB39_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot14 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot14)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()


### OB37 & OB39

# Step 1: Read the data from the CSV file
data <- read.csv("OB37_OB39_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB39") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB37") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB37_OB39_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot15 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot15)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()


### OB38 & OB39

# Step 1: Read the data from the CSV file
data <- read.csv("OB38_OB39_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB39") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB38") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB38-0B39_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot16 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot16)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()

### OB40 & OB42

# Step 1: Read the data from the CSV file
data <- read.csv("OB40_OB42_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB42") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB40") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB40-OB42_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot17 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot17)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()

### OB41 & OB42

# Step 1: Read the data from the CSV file
data <- read.csv("OB41_OB42_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB42") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB41") %>%
  pull(Genus) %>%
  unique()

dev.off()

# Define output file
output_file <- "Fig_OB41-OB42_venn_diagram.jpeg"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1200, res = 150) 

plot18 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera),
  category.names = c(""),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "black"),
  cex = 4,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 1.5,  # Category label size
  cat.fontfamily = "sans"  # Font family for category labels
)

grid.draw(plot18)  # Display the Venn diagram

# Close the graphics device and save the file
dev.off()


### OB35, OB37, OB38, OB39

# Step 1: Read the data from the CSV file
data <- read.csv("OB35_OB37_OB38_OB39_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB39") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB35") %>%
  pull(Genus) %>%
  unique()

sample_type3_genera <- data %>%
  filter(SampleType == "Flukeprint_OB37") %>%
  pull(Genus) %>%
  unique()

sample_type4_genera <- data %>%
  filter(SampleType == "Flukeprint_OB38") %>%
  pull(Genus) %>%
  unique()


# Define output file
output_file <- "Fig_OB35-38_OB39 venn.png"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1000, res = 150)  

# Create the Venn diagram
plot19 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera, 
           SampleType3 = sample_type3_genera, SampleType4 = sample_type4_genera),
  category.names = c("OB39","OB35","OB37","OB38"),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "lightgreen", "lightpink", "violet"),
  cex = 3,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 3,  # Category label size
  cat.fontfamily = "sans",
  margin = 0.1 
)

grid.draw(plot19)  # Draw the Venn diagram


# Close the graphics device and save the file
dev.off()  


### OB23, OB25, OB24

# Step 1: Read the data from the CSV file
data <- read.csv("OB23_OB25_OB24_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB24") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB23") %>%
  pull(Genus) %>%
  unique()

sample_type3_genera <- data %>%
  filter(SampleType == "Flukeprint_OB25") %>%
  pull(Genus) %>%
  unique()


# Define output file
output_file <- "Fig_OB23-25_OB24 venn.png"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1000, res = 150)  

# Create the Venn diagram
plot20 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera, 
           SampleType3 = sample_type3_genera),
  category.names = c("OB24","OB23","OB25"),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "lightgreen", "lightpink"),
  cex = 3,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 3,  # Category label size
  cat.fontfamily = "sans",
  margin = 0.1 
)

grid.draw(plot20)  # Draw the Venn diagram


# Close the graphics device and save the file
dev.off()  


### OB31, OB32, OB33, OB34

# Step 1: Read the data from the CSV file
data <- read.csv("OB31_OB32_OB33_OB34_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB34") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB31") %>%
  pull(Genus) %>%
  unique()

sample_type3_genera <- data %>%
  filter(SampleType == "Flukeprint_OB32") %>%
  pull(Genus) %>%
  unique()

sample_type4_genera <- data %>%
  filter(SampleType == "Flukeprint_OB33") %>%
  pull(Genus) %>%
  unique()

# Define output file
output_file <- "Fig_OB31-33_OB34_venn.png"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1000, res = 150)  

# Create the Venn diagram
plot21 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera, 
           SampleType3 = sample_type3_genera,SampleType4 = sample_type4_genera),
  category.names = c("OB34","OB31","OB32","OB33"),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "lightgreen", "lightpink","violet"),
  cex = 3,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 3,  # Category label size
  cat.fontfamily = "sans",
  margin = 0.1 
)

grid.draw(plot21)  # Draw the Venn diagram


# Close the graphics device and save the file
dev.off()  


### OB40, OB41, OB43

# Step 1: Read the data from the CSV file
data <- read.csv("OB40_OB41_OB42_venn.csv")

# Step 2: Extract unique genera for each sample type
sample_type1_genera <- data %>%
  filter(SampleType == "Seawater Control_OB42") %>%
  pull(Genus) %>%
  unique()

sample_type2_genera <- data %>%
  filter(SampleType == "Flukeprint_OB40") %>%
  pull(Genus) %>%
  unique()

sample_type3_genera <- data %>%
  filter(SampleType == "Flukeprint_OB41") %>%
  pull(Genus) %>%
  unique()


# Define output file
output_file <- "Fig_OB40-41_OB42_venn.png"

# Open a graphics device to save the plot
png(output_file, width = 1400, height = 1000, res = 150)  

# Create the Venn diagram
plot22 <- venn.diagram(
  x = list(SampleType1 = sample_type1_genera, SampleType2 = sample_type2_genera, 
           SampleType3 = sample_type3_genera),
  category.names = c("OB42","OB40","OB41"),
  filename = NULL,  # Output to a variable instead of a file
  output = TRUE,
  # Custom colors: blue for Seawater Control, black for Flukeprint
  fill = c("lightblue", "lightgreen", "lightpink"),
  cex = 3,  # Text size
  fontfamily = "sans",  # Font family
  cat.cex = 3,  # Category label size
  cat.fontfamily = "sans",
  margin = 0.1 
)

grid.draw(plot22)  # Draw the Venn diagram


# Close the graphics device and save the file
dev.off()  


####### END #######