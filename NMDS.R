### OrcaBiome Paper Analysis
### NMDS plots
### March 16, 2025

#Load libraries
library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(goeveg) # scree
library(plyr) # ddply
library(viridis) # viridis colours
library(gridExtra) # arrange plots
library(grid) # arrange plots

# Read in data (1. all samples)
#data is the table produced after metaworks pipeline
A <- read.csv(file="OrcaBiome_results_0.8_cut_genera_removed_noNTC_v3.csv", head=TRUE)


# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID","SampleType", "WhaleID","WhaleSex","Population","Age")

#Sample type = blank or flukeprint, whaleID = ID of whale, sex = sex of whale, population = CA/BC
#edit numbers in [] to correspond with what cols have the above attribute headings

# pivot to make esv matrix (pool across verions)
A.2.esv<-reshape2::dcast(A.1, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.2.esv) <- A.2.esv$SampleName
A.2.esv$SampleName <- NULL

#remove columns with only zeros
esv.notnull<-A.2.esv[,colSums(A.2.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]


#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 18618.15 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to presence-absence matrix
rare.mat[rare.mat>0] <-1

# remove sample OB34 (low reads and outlier in the NMDS)
rare.mat <- rare.mat[rownames(rare.mat) != "OB34_BL_OB31-33BL_NA_BC_NA", ]

# Run NMDS for different k values to check stress levels
stress_values <- numeric()  # Store stress values

for (k in 1:6) {  # Check dimensions from 1 to 6
  nmds <- metaMDS(rare.mat, k = k, trymax = 100, autotransform = FALSE)
  stress_values[k] <- nmds$stress
}

# Save scree plot to PNG
png("Scree.png", width = 800, height = 600, res = 150)
plot(1:6, stress_values, type = "b", pch = 16, col = "blue",
     xlab = "Number of Dimensions (k)", ylab = "Stress",
     main = "NMDS Stress Plot")
abline(h = 0.2, col = "red", lty = 2)  # Threshold for acceptable stress
dev.off()

# Choose K = 2 (stress = 0.06)
nmds2<-metaMDS(rare.mat, k=2, trymax=100)

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n", main="SSU")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.99

# Create grouping matrix for samples by grabbing row names from above matrix
names<-data.frame(row.names(rare.mat), stringsAsFactors = FALSE)

# Rename the column
names(names)<-"sample"

# Copy column to row names
row.names(names)<-names$sample

# Split first column into their own fields
names.1<-data.frame(names, do.call(rbind, strsplit(names$sample,'_')), stringsAsFactors = FALSE)
names(names.1)[2:7]<-c("SampleID","SampleType", "WhaleID","WhaleSex","Population","Age")
#edit numbers in [] to correspond with what cols have the above attribute headings

# Remove first column
names.1 <- names.1[,-1]

# Grab population/species scores from NMDS output
df <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df for ggplot
gg <- merge(df, names.1, by="row.names")

# create factors
gg$SampleType <- factor(gg$SampleType,
                       levels = c("BL", "FP"),
                       labels = c("Seawater Control", "Flukeprint"))
gg$Population <- factor(gg$Population,
                     levels = c("MB", "BC"),
                     labels = c("California", "British Columbia"))
gg$WhaleSex <- factor(gg$WhaleSex,
                      levels = c("M", "F","U","NA"),
                      labels = c("Male","Female","Unknown","NA"))
gg$Age <- factor(gg$Age,
                      levels = c("SA", "CA","OA","NA"),
                      labels = c("Young Adult","Calf","Old Adult","NA"))

# color by population
chulls.pop <- ddply(gg, .(Population), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by population
p1 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.pop, aes(x=NMDS1, y=NMDS2, fill=Population), alpha=0.5) +
  geom_point(data=gg, aes(color=Population)) +
  ggtitle("Bacterial ESVs by Geographic Location") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 20.8, p < 0.01"), 
           hjust = 1, vjust = 1, size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_color_manual(values = c("California" = "#3B4D93", "British Columbia" = "#2A9D8F"),
                     labels = c("California" = "California", "British Columbia" = "British\nColumbia")) +
  scale_fill_manual(values = c("California" = "#3B4D93", "British Columbia" = "#2A9D8F"),
                    labels = c("California" = "California", "British Columbia" = "British\nColumbia"))


p1

#ggsave("Fig_MB vs BC NMDS.jpeg", plot = p1, dpi = 300, height = 6, width = 8, units = "in")

# color by sample type
chulls.sample <- ddply(gg, .(SampleType), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by sample type
p2 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.sample, aes(x=NMDS1, y=NMDS2, fill=SampleType), alpha=0.5) +
  geom_point(data=gg, aes(color=SampleType)) +
  ggtitle("Bacterial ESVs by Sample Type") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 0.90  p = 0.838"), 
           hjust = 1, vjust = 1, size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_color_manual(values = c("Seawater Control" = "lightblue", "Flukeprint" = "black"),
                     labels = c("Seawater Control" = "Seawater\nControl", "Flukeprint" = "Flukeprint")) +
  scale_fill_manual(values = c("Seawater Control" = "lightblue", "Flukeprint" = "black"),
                    labels = c("Seawater Control" = "Seawater\nControl", "Flukeprint" = "Flukeprint"))

p2

#ggsave("Fig_FP vs BL NMDS.jpeg", plot = p2, dpi = 300, height = 6, width = 8, units = "in")

# color by whale sex
chulls.sex <- ddply(gg, .(WhaleSex), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale sex
p3 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.sex, aes(x=NMDS1, y=NMDS2, fill=WhaleSex), alpha=0.5) +
  geom_point(data=gg, aes(color=WhaleSex)) +
  ggtitle("Bacterial ESVs by Whale Sex") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 0.87, p = 0.66"), 
           hjust = 1, vjust = 1, size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_text(size=10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_color_manual(values = c("Male" = "#003C75", "Female" = "#9C2A6C", "Unknown" = "#F04E23", "NA" = "lightgrey" )) +
  scale_fill_manual(values = c("Male" = "#003C75", "Female" = "#9C2A6C", "Unknown" = "#F04E23", "NA" = "lightgrey"))

p3


# color by whale age
chulls.age <- ddply(gg, .(Age), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale age
p4 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.age, aes(x=NMDS1, y=NMDS2, fill=Age), alpha=0.5) +
  geom_point(data=gg, aes(color=Age)) +
  ggtitle("Bacterial ESVs by Whale Age Class") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 1.2, p = 0.135"), 
           hjust = 1, vjust = 1, size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_text(size=10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_color_manual(values = c("Young Adult" = "mediumpurple", "Calf" = "forestgreen", "Old Adult" = "indianred", "NA" = "lightgrey" )) +
  scale_fill_manual(values = c("Young Adult" = "mediumpurple", "Calf" = "forestgreen", "Old Adult" = "indianred", "NA" = "lightgrey"))

p4



plots <- grid.arrange(
  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,   # Arrange plots in a single column
    heights = c(2, 2)  # Adjust heights if necessary
  ),
  left = textGrob("NMDS2", rot = 90, vjust = 1, gp = gpar(fontsize = 12, fontface = "bold")),
  bottom = textGrob("NMDS1", vjust = 0.8, gp = gpar(fontsize = 12, fontface = "bold"))
)


ggsave("Fig4_NMDS.jpeg", plot = plots, dpi = 300, height = 20, width = 26, units = "cm")


### Stats

# Create metadata from rownames 'sample'
env <- gg[,c(1,4:9)]


# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using Bray Curtis (Sorensen) dissimilarity
sor<-vegdist(rare.mat, "bray", binary=TRUE)


# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.sampletype<-betadisper(sor, as.factor(env$SampleType))
bd.population<-betadisper(sor, as.factor(env$Population))
bd.sex<-betadisper(sor, as.factor(env$WhaleSex))
bd.age<-betadisper(sor, as.factor(env$Age))

# check for heterogeneity of beta dispersions within groups BALANCED DESIGN
set.seed(1234)
anova(bd.sampletype) # 
          #Df   Sum Sq  Mean Sq F value Pr(>F)
#Groups     1 0.004314 0.0043139  0.4242 0.5203
#Residuals 27 0.274550 0.0101685  

anova(bd.population) # 
          #Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     1 0.0099127 0.0099127  9.4182 0.004849 **
#Residuals 27 0.0284176 0.0010525 

anova(bd.sex) # 
          #Df   Sum Sq   Mean Sq F value Pr(>F)  
#Groups     3 0.079134 0.0263781  2.9719 0.05098 .
#Residuals 25 0.221893 0.0088757  

anova(bd.age)
          #Df  Sum Sq  Mean Sq F value   Pr(>F)   
#Groups     3 0.19183 0.063943  6.3397 0.002399 **
#Residuals 25 0.25215 0.010086 

pdf("BetaDispersion.pdf")
par(mfrow=c(2,2))
boxplot(bd.sampletype, main="Sample Type")
boxplot(bd.population, main="Population")
boxplot(bd.sex, main="Sex")
boxplot(bd.age, main="Age")
dev.off()

# Use ADONIS to test significance of groupings 
adonis2(sor~SampleType, data=env, permutations=999, strata=env$Population)
        # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Model     1   0.1802 0.03221 0.8985  0.838
#Residual 27   5.4144 0.96779              
#Total    28   5.5946 1.00000

# No significant difference in FP and BL within either pop


adonis2(sor~Population, data=env, permutations=999, strata=env$SampleType)

          #Df SumOfSqs      R2      F Pr(>F)    
#Model     1   2.4355 0.43533 20.816  0.001 ***
#Residual 27   3.1591 0.56467                  
#Total    28   5.5946 1.00000 


# Significant difference in California vs BC (by sample type)


adonis2(sor~WhaleSex, data=env, permutations=999, strata=env$Population)

          #Df SumOfSqs      R2      F Pr(>F)
#Model     3   0.5275 0.0943 0.8676   0.66
#Residual 25   5.0670 0.9057              
#Total    28   5.5946 1.0000 

# No significant difference in whale sex within either pop


adonis2(sor~Age, data=env, permutations=999, strata=env$Population)

          #Df SumOfSqs      R2      F Pr(>F)
#Model     3   0.7187 0.12846 1.2283  0.135
#Residual 25   4.8759 0.87154              
#Total    28   5.5946 1.00000   


# Read in just flukeprint bacteria (both populations), test for sex and age)


# Read the file containing the flukeprints families
flukeprints_ESV_file <- "ESV_FP.csv"  # Replace with your actual file path
flukeprints_ESV <- read.csv(flukeprints_ESV_file)$GlobalESV


# Filter long_data to only include the ESVs found in flukeprints
FP_ESV <- dplyr::filter(A, GlobalESV %in% flukeprints_ESV)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(FP_ESV, do.call(rbind, str_split(FP_ESV$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID","SampleType", "WhaleID","WhaleSex","Population","Age")

#Sample type = blank or flukeprint, whaleID = ID of whale, sex = sex of whale, population = CA/BC
#edit numbers in [] to correspond with what cols have the above attribute headings

# pivot to make esv matrix (pool across verions)
A.2.esv<-reshape2::dcast(A.1, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.2.esv) <- A.2.esv$SampleName
A.2.esv$SampleName <- NULL

#remove columns with only zeros
esv.notnull<-A.2.esv[,colSums(A.2.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]


#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 130

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to presence-absence matrix
rare.mat[rare.mat>0] <-1

# Run NMDS for different k values to check stress levels
stress_values <- numeric()  # Store stress values

for (k in 1:6) {  # Check dimensions from 1 to 6
  nmds <- metaMDS(rare.mat, k = k, trymax = 100, autotransform = FALSE)
  stress_values[k] <- nmds$stress
}

# Save scree plot to PNG
png("Scree_FP ESVs only.png", width = 800, height = 600, res = 150)
plot(1:6, stress_values, type = "b", pch = 16, col = "blue",
     xlab = "Number of Dimensions (k)", ylab = "Stress",
     main = "NMDS Stress Plot")
abline(h = 0.2, col = "red", lty = 2)  # Threshold for acceptable stress
dev.off()

# Choose K = 2 (stress = 0.06)
nmds2<-metaMDS(rare.mat, k=2, trymax=100)

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_FP ESVs only.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n", main="SSU")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.99

# Create grouping matrix for samples by grabbing row names from above matrix
names<-data.frame(row.names(rare.mat), stringsAsFactors = FALSE)

# Rename the column
names(names)<-"sample"

# Copy column to row names
row.names(names)<-names$sample

# Split first column into their own fields
names.1<-data.frame(names, do.call(rbind, strsplit(names$sample,'_')), stringsAsFactors = FALSE)
names(names.1)[2:7]<-c("SampleID","SampleType", "WhaleID","WhaleSex","Population","Age")
#edit numbers in [] to correspond with what cols have the above attribute headings

# Remove first column
names.1 <- names.1[,-1]

# Grab population/species scores from NMDS output
df <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df for ggplot
gg <- merge(df, names.1, by="row.names")

# create factors
gg$Population <- factor(gg$Population,
                        levels = c("MB", "BC"),
                        labels = c("California", "British Columbia"))
gg$WhaleSex <- factor(gg$WhaleSex,
                      levels = c("M", "F","U"),
                      labels = c("Male","Female","Unknown"))
gg$Age <- factor(gg$Age,
                 levels = c("SA","CA","OA","NA"),
                 labels = c("Young Adult","Calf","Old Adult","NA"))


# color by population
chulls.pop <- ddply(gg, .(Population), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by population
p9 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.pop, aes(x=NMDS1, y=NMDS2, fill=Population), alpha=0.5) +
  geom_point(data=gg, aes(color=Population)) +
  ggtitle("Unique Flukeprint Bacterial ESVs by Geographic Location") +
  annotate("text", x = max(gg$NMDS1), y = min(gg$NMDS2), 
           label = paste0("F = 2.1, p < 0.01 "), 
           hjust = 0.9, vjust = 0.4, size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_color_manual(values = c("California" = "#3B4D93", "British Columbia" = "#2A9D8F"),
                     labels = c("California" = "California", "British Columbia" = "British\nColumbia")) +
  scale_fill_manual(values = c("California" = "#3B4D93", "British Columbia" = "#2A9D8F"),
                    labels = c("California" = "California", "British Columbia" = "British\nColumbia"))


p9

#ggsave("Fig_MB vs BC NMDS.jpeg", plot = p1, dpi = 300, height = 6, width = 8, units = "in")


# color by whale sex
chulls.sex <- ddply(gg, .(WhaleSex), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale sex
p10 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.sex, aes(x=NMDS1, y=NMDS2, fill=WhaleSex), alpha=0.5) +
  geom_point(data=gg, aes(color=WhaleSex)) +
  ggtitle("Unique Flukeprint Bacterial ESVs by Whale Sex") +
  annotate("text", x = max(gg$NMDS1), y = min(gg$NMDS2), 
           label = paste0("F = 1.1, p = 0.105 "), 
           hjust = 0.9, vjust = 0.4, size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_color_manual(values = c("Male" = "#003C75", "Female" = "#9C2A6C", "Unknown" = "#F04E23", "NA" = "lightgrey" )) +
  scale_fill_manual(values = c("Male" = "#003C75", "Female" = "#9C2A6C", "Unknown" = "#F04E23", "NA" = "lightgrey"))

p10


# color by whale age
chulls.age <- ddply(gg, .(Age), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale age
p11 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.age, aes(x=NMDS1, y=NMDS2, fill=Age), alpha=0.5) +
  geom_point(data=gg, aes(color=Age)) +
  ggtitle("Unique Flukeprint Bacterial ESVs by Whale Age Class") +
  annotate("text", x = max(gg$NMDS1), y = min(gg$NMDS2), 
           label = paste0("F = 1.0, p = 0.451 "), 
           hjust = 0.9, vjust = 0.4, size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_text(size=10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_color_manual(values = c("Young Adult" = "mediumpurple", "Calf" = "forestgreen", "Old Adult" = "indianred", "NA" = "lightgrey" )) +
  scale_fill_manual(values = c("Young Adult" = "mediumpurple", "Calf" = "forestgreen", "Old Adult" = "indianred", "NA" = "lightgrey"))

p11

plots <- grid.arrange(
  arrangeGrob(
    p9, p10, p11,
    ncol = 1,   # Arrange plots in a single column
    heights = c(2,2,2)  # Adjust heights if necessary
  ),
  left = textGrob("NMDS2", rot = 90, vjust = 1, gp = gpar(fontsize = 12, fontface = "bold")),
  bottom = textGrob("NMDS1", vjust = 0.8, gp = gpar(fontsize = 12, fontface = "bold"))
)

ggsave("Fig5_NMDS_FP only.jpeg", plot = plots, dpi = 300, height = 24, width = 18, units = "cm")

### Stats

# Create metadata from rownames 'sample'
env <- gg[,c(1,4:9)]


# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using Bray Curtis (Sorensen) dissimilarity
sor<-vegdist(rare.mat, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.population<-betadisper(sor, as.factor(env$Population))
bd.sex<-betadisper(sor, as.factor(env$WhaleSex))
bd.age<-betadisper(sor, as.factor(env$Age))

# check for heterogeneity of beta dispersions within groups BALANCED DESIGN
set.seed(1234)

anova(bd.population) # 
#Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     1 0.0008486 0.00084864  0.5023 0.4887
#Residuals 16 0.0270325 0.00168953    

anova(bd.sex) # 
#Df   Sum Sq   Mean Sq F value Pr(>F)  
#Groups     2 0.044689 0.0223445  14.059 0.0003636 ***
#Residuals 15 0.023839 0.0015893     

anova(bd.age)
          #Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     2 0.37619 0.188093  109.28 1.142e-09 ***
#Residuals 15 0.02582 0.001721 

pdf("BetaDispersion_unique FP only.pdf")
par(mfrow=c(2,2))
boxplot(bd.sampletype, main="Sample Type")
boxplot(bd.population, main="Population")
boxplot(bd.sex, main="Sex")
boxplot(bd.age, main="Age")
dev.off()

adonis2(sor~Population, data=env, permutations=999)

          #Df SumOfSqs      R2      F Pr(>F)    
#Model     1   0.8742 0.11372 2.053  0.001 ***
#Residual 16   6.8129 0.88628                 
#Total    17   7.6870 1.00000   


# Significant difference in California vs BC


adonis2(sor~WhaleSex, data=env, permutations=999)

        #Df SumOfSqs      R2      F Pr(>F)
#Model     2   1.0175 0.13237 1.1442  0.105
#Residual 15   6.6695 0.86763              
#Total    17   7.6870 1.00000 

# No significant difference in whale sex within either pop

adonis2(sor~Age, data=env, permutations=999)

        #Df SumOfSqs      R2      F Pr(>F)  
#Model     2   0.9087 0.11821 1.0054  0.451
#Residual 15   6.7784 0.88179              
#Total    17   7.6870 1.00000 

# No significant difference in whale age within either pop

############## END ##################







