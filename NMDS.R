### OrcaBiome Paper Analysis
### NMDS plots
### December 02, 2025

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

# Read in data
A <- read.csv(file="OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", head=TRUE)

# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID","Type", "WhaleID","WhaleSex","Population","Age")

# Pivot to make esv matrix (pool across verions)
A.2.esv<-reshape2::dcast(A.1, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Move sample to rownames then delete
rownames(A.2.esv) <- A.2.esv$SampleName
A.2.esv$SampleName <- NULL

# Remove columns with only zeros
esv.notnull<-A.2.esv[,colSums(A.2.esv) !=0]

# Remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# Calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 9892.75 

# Set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# Remove sample OB34 (low reads and outlier in the NMDS)
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
names(names.1)[2:7]<-c("SampleID","Type", "WhaleID","WhaleSex","Population","Age")
#edit numbers in [] to correspond with what cols have the above attribute headings

# Remove first column
names.1 <- names.1[,-1]

# Grab population/species scores from NMDS output
df <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df for ggplot
gg <- merge(df, names.1, by="row.names")

# Create factors
gg$Type <- factor(gg$Type,
                       levels = c("BL", "FP"),
                       labels = c("Seawater Control", "Flukeprint"))
gg$Population <- factor(gg$Population,
                     levels = c("MB", "BC"),
                     labels = c("California", "British Columbia"))
gg$WhaleSex <- factor(gg$WhaleSex,
                      levels = c("M","F","U","NA"),
                      labels = c("Male","Female","Unknown","NA"))
gg$Age <- factor(gg$Age,
                      levels = c("SA","CA","OA","NA"),
                      labels = c("Young Adult","Calf","Old Adult","NA"))

# color by population
chulls.pop <- ddply(gg, .(Population), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by population
p1 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.pop, aes(x=NMDS1, y=NMDS2, fill=Population), alpha=0.5) +
  geom_point(data=gg, aes(color=Population)) +
  ggtitle("Bacterial ESVs by Geographic Location") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 23.8, p < 0.001"), 
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

# Color by sample type
chulls.sample <- ddply(gg, .(Type), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by sample type
p2 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.sample, aes(x=NMDS1, y=NMDS2, fill=Type), alpha=0.5) +
  geom_point(data=gg, aes(color=Type)) +
  ggtitle("Bacterial ESVs by Sample Type") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 0.92  p = 0.761"), 
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

# Color by whale sex
chulls.sex <- ddply(gg, .(WhaleSex), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale sex
p3 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.sex, aes(x=NMDS1, y=NMDS2, fill=WhaleSex), alpha=0.5) +
  geom_point(data=gg, aes(color=WhaleSex)) +
  ggtitle("Bacterial ESVs by Whale Sex") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 0.89, p = 0.406"), 
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


# Color by whale age
chulls.age <- ddply(gg, .(Age), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale age
p4 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.age, aes(x=NMDS1, y=NMDS2, fill=Age), alpha=0.5) +
  geom_point(data=gg, aes(color=Age)) +
  ggtitle("Bacterial ESVs by Whale Age Class") +
  annotate("text", x = max(gg$NMDS1), y = max(gg$NMDS2), 
           label = paste0("F = 1.31, p = 0.192"), 
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

plots

ggsave("FigS11_NMDS.jpeg", plot = plots, dpi = 300, height = 20, width = 26, units = "cm")


### Stats

# Create metadata from rownames 'sample'
env <- gg[,c(1,4:9)]

# Assess dispersion (variance) using ANOVA
# Create distance matrix 
sor<-vegdist(rare.mat, "bray", binary=FALSE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.type <- betadisper(sor, as.factor(env$Type))
bd.population <- betadisper(sor, as.factor(env$Population))
bd.sex <- betadisper(sor, as.factor(env$WhaleSex))
bd.age <- betadisper(sor, as.factor(env$Age))

# check for heterogeneity of beta dispersions within groups BALANCED DESIGN
set.seed(1234)
anova(bd.type) # 
          #Df   Sum Sq  Mean Sq F value Pr(>F)
#Groups     1 0.03327 0.033270  0.9665 0.3343
#Residuals 27 0.92939 0.034422    

anova(bd.population) # 
          #Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     1 0.14593 0.145929  11.378 0.00226 **
#Residuals 27 0.34630 0.012826   

anova(bd.sex) # 
          #Df   Sum Sq   Mean Sq F value Pr(>F)  
#Groups     3 0.19233 0.064110  2.1803 0.1155
#Residuals 25 0.73511 0.029405  

anova(bd.age)
          #Df  Sum Sq  Mean Sq F value   Pr(>F)   
#Groups     3 0.29217 0.097389  2.8438 0.05805 .
#Residuals 25 0.85616 0.034246   

pdf("BetaDispersion.pdf")
par(mfrow=c(2,2))
boxplot(bd.type, main="Sample Type")
boxplot(bd.population, main="Population")
boxplot(bd.sex, main="Sex")
boxplot(bd.age, main="Age")
dev.off()

# Use ADONIS to test significance of groupings 
adonis2(sor~Type, data=env, permutations=999, strata=env$Population)
        # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Model     1   0.2328 0.03288 0.918  0.761
#Residual 27   6.8474 0.96712             
#Total    28   7.0802 1.00000  

# No significant difference in FP and BL within either pop


adonis2(sor~Population, data=env, permutations=999, strata=env$Type)

          #Df SumOfSqs      R2      F Pr(>F)    
#Model     1   3.3168 0.46846 23.796  0.001 ***
#Residual 27   3.7634 0.53154                  
#Total    28   7.0802 1.00000   


# Significant difference in California vs BC (by sample type)


adonis2(sor~WhaleSex, data=env, permutations=999, strata=env$Population)

          #Df SumOfSqs      R2      F Pr(>F)
#Model     3   0.6839 0.09659 0.891  0.406
#Residual 25   6.3963 0.90341             
#Total    28   7.0802 1.00000   

# No significant difference in whale sex within either pop


adonis2(sor~Age, data=env, permutations=999, strata=env$Population)

          #Df SumOfSqs      R2      F Pr(>F)
#Model     3   0.9642 0.13618 1.3137  0.192
#Residual 25   6.1160 0.86382              
#Total    28   7.0802 1.00000      


####### Read in just flukeprint bacteria (both populations), test for sex and age)


# Read the file containing the flukeprints families
flukeprints_ESV_file <- "ESV_FP.csv"  # Replace with your actual file path
flukeprints_ESV <- read.csv(flukeprints_ESV_file)$GlobalESV


# Filter long_data to only include the ESVs found in flukeprints
FP_ESV <- dplyr::filter(A, GlobalESV %in% flukeprints_ESV)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(FP_ESV, do.call(rbind, str_split(FP_ESV$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[26:31] <- c("SampleID","Type", "WhaleID","WhaleSex","Population","Age")

# Pivot to make esv matrix (pool across verions)
A.2.esv<-reshape2::dcast(A.1, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Move sample to rownames then delete
rownames(A.2.esv) <- A.2.esv$SampleName
A.2.esv$SampleName <- NULL

#Remove columns with only zeros
esv.notnull<-A.2.esv[,colSums(A.2.esv) !=0]

#Remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]


#Calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 127

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

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
names(names.1)[2:7]<-c("SampleID","Type", "WhaleID","WhaleSex","Population","Age")
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


# Color by population
chulls.pop <- ddply(gg, .(Population), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by population
p9 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.pop, aes(x=NMDS1, y=NMDS2, fill=Population), alpha=0.5) +
  geom_point(data=gg, aes(color=Population)) +
  ggtitle("Unique Flukeprint Bacterial ESVs by Geographic Location") +
  annotate("text", x = max(gg$NMDS1), y = min(gg$NMDS2), 
           label = paste0("F = 1.67, p < 0.001 "), 
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


# Color by whale sex
chulls.sex <- ddply(gg, .(WhaleSex), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale sex
p10 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.sex, aes(x=NMDS1, y=NMDS2, fill=WhaleSex), alpha=0.5) +
  geom_point(data=gg, aes(color=WhaleSex)) +
  ggtitle("Unique Flukeprint Bacterial ESVs by Whale Sex") +
  annotate("text", x = max(gg$NMDS1), y = min(gg$NMDS2), 
           label = paste0("F = 1.10, p = 0.095 "), 
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


# Color by whale age
chulls.age <- ddply(gg, .(Age), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by whale age
p11 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.age, aes(x=NMDS1, y=NMDS2, fill=Age), alpha=0.5) +
  geom_point(data=gg, aes(color=Age)) +
  ggtitle("Unique Flukeprint Bacterial ESVs by Whale Age Class") +
  annotate("text", x = max(gg$NMDS1), y = min(gg$NMDS2), 
           label = paste0("F = 1.01, p = 0.381 "), 
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

ggsave("Fig4_NMDS_FP only.jpeg", plot = plots, dpi = 300, height = 24, width = 18, units = "cm")

### Stats

# Create metadata from rownames 'sample'
env <- gg[,c(1,4:9)]


# Assess dispersion (variance) using ANOVA
# Create distance matrix
sor<-vegdist(rare.mat, "bray", binary=FALSE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.population<-betadisper(sor, as.factor(env$Population))
bd.sex<-betadisper(sor, as.factor(env$WhaleSex))
bd.age<-betadisper(sor, as.factor(env$Age))

# check for heterogeneity of beta dispersions within groups BALANCED DESIGN
set.seed(1234)

anova(bd.population) # 
#Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     1 0.0013348 0.00133479  1.5005 0.2383
#Residuals 16 0.0142326 0.00088954    

anova(bd.sex) # 
#Df   Sum Sq   Mean Sq F value Pr(>F)  
#Groups     2 0.044689 0.0223445  14.059 0.0003636 ***
#Residuals 15 0.023839 0.0015893     

anova(bd.age)
          #Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     2 0.050086 0.0250429  67.887 3.043e-08 ***
#Residuals 15 0.005533 0.0003689    

pdf("BetaDispersion_unique FP only.pdf")
par(mfrow=c(2,2))
boxplot(bd.sampletype, main="Sample Type")
boxplot(bd.population, main="Population")
boxplot(bd.sex, main="Sex")
boxplot(bd.age, main="Age")
dev.off()

adonis2(sor~Population, data=env, permutations=999)

          #Df SumOfSqs      R2      F Pr(>F)    
#Model     1   0.7578 0.09455 1.6708  0.001 ***
#Residual 16   7.2570 0.90545                  
#Total    17   8.0148 1.00000 


# Significant difference in California vs BC


adonis2(sor~WhaleSex, data=env, permutations=999)

        #Df SumOfSqs      R2      F Pr(>F)
#Model     2   1.0233 0.12768 1.0978  0.095 .
#Residual 15   6.9915 0.87232                
#Total    17   8.0148 1.00000   

# No significant difference in whale sex within either pop

adonis2(sor~Age, data=env, permutations=999)

        #Df SumOfSqs      R2      F Pr(>F)  
#Model     2   0.9537 0.11899 1.0129  0.381
#Residual 15   7.0611 0.88101              
#Total    17   8.0148 1.00000  

# No significant difference in whale age within either pop

############## END ##################







