### Orcabiome eDNA-Skin Models ###
#Created by: Brittany Visona-Kelly
#Created: June 6, 2025
#Edited: December 04, 2025

###########################
## Data Processing
###########################

# Load libaries
library(DescTools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(tidyr)
library(performance)
library(effects)
library(ggeffects)

biome <- read.csv(("Orcabiome data eDNA models.csv"), stringsAsFactors =TRUE, na.strings = c(""))
str(biome)

biome$Age <- 2024 - as.numeric(as.character(biome$Birth.year))
biome$Age

toJulianDate <- function(dates) {
  dates <- as.Date(dates)
  
  sapply(dates, function(date) {
    y <- as.numeric(format(date, "%Y"))
    m <- as.numeric(format(date, "%m"))
    d <- as.numeric(format(date, "%d"))
    
    if (m <= 2) {
      y <- y - 1
      m <- m + 12
    }
    
    A <- floor(y / 100)
    B <- 2 - A + floor(A / 4)
    
    JD <- floor(365.25 * (y + 4716)) +
      floor(30.6001 * (m + 1)) +
      d + B - 1524.5
    
    return(JD)
  })
}

biome$Sample.dateformat <- as.Date(as.character(biome$Sample.date), format = "%m/%d/%Y")
biome$Sample.juliandate<- toJulianDate(biome$Sample.dateformat)
biome$Sample.juliandate

biome$Month<-as.factor(format(biome$Sample.dateformat, "%B"))

biome$Image.dateformat <- as.Date(as.character(biome$Image.date), format = "%d-%b-%y")
biome$Image.juliandate<- toJulianDate(biome$Image.dateformat)
biome$Image.juliandate

biome$diff.sampleimagedate <- biome$Sample.juliandate-biome$Image.juliandate
biome$diff.sampleimagedate

mean(abs(biome$diff.sampleimagedate))
#[1] 21.46667
range(abs(biome$diff.sampleimagedate))
#[1]   0 116


Mode(abs(biome$diff.sampleimagedate))
#[1] 0
#attr(,"freq")
#[1] 9

#Remove section 1 images with a score =1
biome$S1.score<-as.factor(biome$S1.score)
biome$S1.lesions.pixels.avg<-as.numeric(as.character(biome$S1.lesions.pixels.avg))

biome$S1.lesions.pixels.avg.omit1 <- ifelse(
  biome$S1.score == "1",
  NA,
  biome$S1.lesions.pixels.avg
)
biome$S1.lesions.pixels.avg.omit1 

biome$S1.allpixels.avg.omit1 <- ifelse(
  biome$S1.score == "1" & !is.na(biome$S1.score),
  NA,
  biome$S1.allpixels.avg 
)
biome$S1.allpixels.avg.omit1 #4 NAs created (4 images with S1.score =1)

biome[c(biome$S1.score == "1"), ] #Check number of images where S1.score = 1

###########################
## Average per whale
###########################

str(biome)

subset_biome <- biome %>%
  select(Whale.ID, Region, Age, Sex, Month, Number.likely.lesion.types, diff.sampleimagedate, S1.lesions.pixels.avg.omit1, S2.lesions.pixels.avg, 
         S3.lesions.pixels.avg, S4.lesions.pixels.avg, S1.allpixels.avg.omit1,S2.allpixels.avg, S3.allpixels.avg, S4.allpixels.avg, 
         S1.scar.score.avg, S2.scar.score.avg, S3.scar.score.avg, S4.scar.score.avg)
str(subset_biome)

subset_biome <- subset_biome %>%
  mutate(across(
    c(S2.lesions.pixels.avg, S3.lesions.pixels.avg, S4.lesions.pixels.avg,
      S1.allpixels.avg.omit1, S2.allpixels.avg, S3.allpixels.avg, S4.allpixels.avg,
      S1.scar.score.avg, S2.scar.score.avg, S3.scar.score.avg, S4.scar.score.avg),
    ~ suppressWarnings(as.numeric(as.character(.)))
  ))
str(subset_biome)

biome.whale <- subset_biome %>% group_by(Whale.ID) %>%
  summarise(
    S1.lesions.pixels.avg.omit1 = sum(na.omit(S1.lesions.pixels.avg.omit1)),
    S2.lesions.pixels.avg = sum(na.omit(S2.lesions.pixels.avg)),
    S3.lesions.pixels.avg = sum(na.omit(S3.lesions.pixels.avg)),
    S4.lesions.pixels.avg  = sum(na.omit(S4.lesions.pixels.avg)),
    S1.allpixels.avg.omit1  = sum(na.omit(S1.allpixels.avg.omit1)),
    S2.allpixels.avg = sum(na.omit(S2.allpixels.avg)),
    S3.allpixels.avg = sum(na.omit(S3.allpixels.avg)),
    S4.allpixels.avg  = sum(na.omit(S4.allpixels.avg)),
    S1.scar.score.avg  = mean(na.omit(S1.scar.score.avg)),
    S2.scar.score.avg  = mean(na.omit(S2.scar.score.avg)),
    S3.scar.score.avg  = mean(na.omit(S3.scar.score.avg)),
    S4.scar.score.avg  = mean(na.omit(S4.scar.score.avg)), 
  )
str(biome.whale)
biome.whale

#Total lesion percentage per whale
biome.whale$lesions.pixels.total <- rowSums(
  cbind(
    biome.whale$S1.lesions.pixels.avg.omit1,
    biome.whale$S2.lesions.pixels.avg,
    biome.whale$S3.lesions.pixels.avg,
    biome.whale$S4.lesions.pixels.avg),
  na.rm = TRUE
)
biome.whale$all.pixels.total <- rowSums(
  cbind(
    biome.whale$S1.allpixels.avg.omit1,
    biome.whale$S2.allpixels.avg,
    biome.whale$S3.allpixels.avg,
    biome.whale$S4.allpixels.avg),
  na.rm = TRUE
)

biome.whale$percentlesions.total <- (biome.whale$lesions.pixels.total / biome.whale$all.pixels.total)*100
biome.whale$percentlesions.total

#Convert mean scar score to a factor

for (i in 1:nrow(biome.whale)) {
  score <- biome.whale$S1.scar.score.avg[i]
  if (!is.na(score)) {
    if (score == 1) {
      biome.whale$S1.scarscore[i] <- "None"
    } else if (score > 1 & score <= 2.25) {
      biome.whale$S1.scarscore[i] <- "Mild"
    } else if (score > 2.25 & score <= 3.25) {
      biome.whale$S1.scarscore[i] <- "Moderate"
    } else if (score > 3.25) {
      biome.whale$S1.scarscore[i] <- "Severe"
    }
  } else {
    biome.whale$S1.scarscore[i] <- NA  # optional: set NA explicitly
  }
}
biome.whale$S1.scarscore

for (i in 1:nrow(biome.whale)) {
  score <- biome.whale$S2.scar.score.avg[i]
  if (!is.na(score)) {
    if (score == 1) {
      biome.whale$S2.scarscore[i] <- "None"
    } else if (score > 1 & score <= 2.25) {
      biome.whale$S2.scarscore[i] <- "Mild"
    } else if (score > 2.25 & score <= 3.25) {
      biome.whale$S2.scarscore[i] <- "Moderate"
    } else if (score > 3.25) {
      biome.whale$S2.scarscore[i] <- "Severe"
    }
  } else {
    biome.whale$S2.scarscore[i] <- NA  # optional: set NA explicitly
  }
}
biome.whale$S2.scarscore

for (i in 1:nrow(biome.whale)) {
  score <- biome.whale$S3.scar.score.avg[i]
  if (!is.na(score)) {
    if (score == 1) {
      biome.whale$S3.scarscore[i] <- "None"
    } else if (score > 1 & score <= 2.25) {
      biome.whale$S3.scarscore[i] <- "Mild"
    } else if (score > 2.25 & score <= 3.25) {
      biome.whale$S3.scarscore[i] <- "Moderate"
    } else if (score > 3.25) {
      biome.whale$S3.scarscore[i] <- "Severe"
    }
  } else {
    biome.whale$S3.scarscore[i] <- NA  # optional: set NA explicitly
  }
}
biome.whale$S3.scarscore

for (i in 1:nrow(biome.whale)) {
  score <- biome.whale$S4.scar.score.avg[i]
  if (!is.na(score)) {
    if (score == 1) {
      biome.whale$S4.scarscore[i] <- "None"
    } else if (score > 1 & score <= 2.25) {
      biome.whale$S4.scarscore[i] <- "Mild"
    } else if (score > 2.25 & score <= 3.25) {
      biome.whale$S4.scarscore[i] <- "Moderate"
    } else if (score > 3.25) {
      biome.whale$S4.scarscore[i] <- "Severe"
    }
  } else {
    biome.whale$S4.scarscore[i] <- NA  # optional: set NA explicitly
  }
}
biome.whale$S4.scarscore

str(biome.whale)

write.csv(biome.whale, "biome.whale.csv") 

### Plot lesion % and scaring by region ###
regions <- read.csv(("regions.csv"), stringsAsFactors =TRUE, na.strings = c(""))
str(regions)
regions$Region <- as.factor(regions$Region)
regions$Location <- factor(regions$Location, levels = c("CA", "BC"))
regions$percentlesions <- (regions$lesions.pixels.avg / regions$allpixels.avg)*100
regions$percentlesions


region.plot <- ggplot(regions, aes(x = Whale.ID, y = percentlesions, fill = Region)) +
  geom_col(position = "stack") +
  labs(y = "Percentage of VSA with Lesions",x = "Whale ID",fill = "Region") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "right") +
  scale_fill_manual(values=c( "#fde725", "#4ac16d", "#365c8d",  "#481b6d","#d0e11c","#a0da39","#73d056","#2db27d", "#1fa187", "#21918c", "#277f8e", 
                              "#2e6e8e","#3f4788", "#46327e",  "#440154")) 
region.plot

regions_label <- regions %>%
  mutate(Region = factor(Region, levels = rev(levels(factor(Region))))) %>%
  filter(scarscore.cat != "None") %>%  # remove 'None' from labels
  group_by(Whale.ID) %>%
  arrange(Region, .by_group = TRUE) %>%
  mutate(
    total = sum(scarscore.freq.freq),
    ymin = cumsum(scarscore.freq.freq) - scarscore.freq.freq,
    ymax = cumsum(scarscore.freq.freq),
    ymid = ymin + scarscore.freq.freq / 2
  )

regionscar.plot <- ggplot(regions, aes(x = Whale.ID, y = scarscore.freq.freq, fill = Region)) +
  geom_col(position = "stack") +
  geom_text(data = regions_label,aes(x = Whale.ID, y = ymid, label = scarscore.cat),color = "white", size = 3) +
  labs(y = "Scar Index of VSA (Relative Frequency)", x = "Whale ID", fill = "Region") +
  theme_bw() +theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position = "right") +
  scale_fill_manual(values = c(
    "#fde725", "#4ac16d", "#365c8d", "#481b6d", "#d0e11c", "#a0da39",
    "#73d056", "#2db27d", "#1fa187", "#21918c", "#277f8e", "#2e6e8e",
    "#3f4788", "#46327e", "#440154"))

regionscar.plot

##### Plot both ###
regionswhale_plot<- ggarrange(region.plot, regionscar.plot, nrow = 2,ncol = 1, 
                              labels = c("(a)", "(b)"), font.label=list(color="black",size=11),
                              common.legend = TRUE, heights = c(4, 4))
regionswhale_plot
#ggexport(regionswhale_plot, filename = "Fig7_regionswhale_plot.pdf", width = 10, height = 10, units = "in", dpi = 800)

###########################
## PCA
###########################

#Combine scare scores for sections 1-4 to reduce number of variables
scar_data <- biome.whale[ , c("S1.scar.score.avg", "S2.scar.score.avg", "S3.scar.score.avg", "S4.scar.score.avg")]
str(scar_data)
scar_data <- data.frame(lapply(scar_data, function(x) as.numeric(as.character(x))))
scar_data[is.na(scar_data)] <- 0
pca_result <- prcomp(scar_data, center = TRUE, scale. = TRUE) # Run PCA (scale = TRUE standardizes variables)
summary(pca_result) 
#Importance of components:
#  PC1    PC2    PC3     PC4
#Standard deviation     1.3781 1.1087 0.7661 0.53370
#Proportion of Variance 0.4748 0.3073 0.1467 0.07121
#Cumulative Proportion  0.4748 0.7821 0.9288 1.00000 #PC1 and PC2 explain 78% of variance

pca_result$rotation
#                         PC1        PC2        PC3        PC4
#S1.scar.score.avg -0.4655342 -0.5206599  0.6015875 -0.3876644
#S2.scar.score.avg -0.4038746 -0.5887116 -0.6421697  0.2791453
#S3.scar.score.avg -0.5287468  0.4856127 -0.3691903 -0.5901743
#S4.scar.score.avg -0.5836010  0.3827692  0.2990137  0.6507599

#Add PCA results to data set
biome.whale$scar_PC1 <- pca_result$x[,1] #Negative PC1 = more scarring across all 4 regions (Scarring amount)
biome.whale$scar_PC2 <- pca_result$x[,2] #Negative PC2 = more scarring S1 and S2, Positive PC2 = more scarring S3 and S4 (Scarring distribution)

###########################
## Data Exploration
###########################
#Add other variables to biome.whale dataframe
biome.whale$Total.ESVs <- biome[match(biome.whale$Whale.ID, biome$Whale.ID),"Total.ESVs", ] 
biome.whale$Total.bacteria.genera <- biome[match(biome.whale$Whale.ID, biome$Whale.ID),"Total.bacteria.genera", ] 
biome.whale$Region<- biome[match(biome.whale$Whale.ID, biome$Whale.ID),"Region", ] 
biome.whale$Age <- biome[match(biome.whale$Whale.ID, biome$Whale.ID),"Age", ] 
biome.whale$Sex <- biome[match(biome.whale$Whale.ID, biome$Whale.ID),"Sex", ] 
biome.whale$Month <- biome[match(biome.whale$Whale.ID, biome$Whale.ID),"Month", ] 
biome.whale$Number.likely.lesion.types <- biome[match(biome.whale$Whale.ID, biome$Whale.ID),"Number.likely.lesion.types", ] 

#Check for collinearity in exp. variables

GGally::ggpairs(biome.whale[, c("Age", "Number.likely.lesion.types","percentlesions.total", "scar_PC1", "scar_PC2")])
#Number.likely.lesion.types and percentlesions.total are correlated (> 0.80). Only use 1 in models. 

#Test normality

ggdensity(biome.whale$Age,xlab = "Age")
ggqqplot(biome.whale$Age,title = "Q-Q Plot of Age")
shapiro.test(biome.whale$Age) #W = 0.83043, p-value = 0.009285 *not normally distributed
ggdensity(biome.whale$Number.likely.lesion.types,xlab = "Number.likely.lesion.types")
ggqqplot(biome.whale$Number.likely.lesion.types,title = "Q-Q Plot of Number.likely.lesion.types")
shapiro.test(biome.whale$Number.likely.lesion.types) #W = 0.86739, p-value = 0.0309 *not normally distributed
ggdensity(biome.whale$percentlesions.total,xlab = "percentlesions.total")
ggqqplot(biome.whale$percentlesions.total,title = "Q-Q Plot of percentlesions.total")
shapiro.test(biome.whale$percentlesions.total)#W = 0.75534, p-value = 0.001038 *not normally distributed
ggdensity(biome.whale$scar_PC1,xlab = "scar_PC1")
ggqqplot(biome.whale$scar_PC1,title = "Q-Q Plot of scar_PC1")
shapiro.test(biome.whale$scar_PC1) #W = 0.94646, p-value = 0.4706
ggdensity(biome.whale$scar_PC2,xlab = "scar_PC2")
ggqqplot(biome.whale$scar_PC2,title = "Q-Q Plot of scar_PC2")
shapiro.test(biome.whale$scar_PC2) #W = 0.96279, p-value = 0.7408
ggdensity(biome.whale$Total.ESVs,xlab = "Total.ESVs")
ggqqplot(biome.whale$Total.ESVs,title = "Q-Q Plot of Total.ESVs")
shapiro.test(biome.whale$Total.ESVs) #W = 0.77074, p-value = 0.001587 *not normally distributed
ggdensity(biome.whale$Total.bacteria.genera,xlab = "Total.bacteria.genera")
ggqqplot(biome.whale$Total.bacteria.genera,title = "Q-Q Plot of Total.bacteria.genera")
shapiro.test(biome.whale$Total.bacteria.genera) #W = 0.91911, p-value = 0.1867

#################################################
## Statistical Analyses - Comparisons
################################################
### Lesions by region ###
biome.whale %>% count(Region, sort = TRUE)
#Region     n
#  1 BC        11
#  2 CA         4

#### Run permutation on small, non-normal dataset ###
# Extract the data
lesion_data <- biome.whale$percentlesions.total
region_group <- biome.whale$Region

# Ensure it's a two-group comparison
if (length(unique(region_group)) != 2) {
  stop("This permutation test only works for two groups.")
}

# Split the data
group1 <- lesion_data[region_group == unique(region_group)[1]]
group2 <- lesion_data[region_group == unique(region_group)[2]]

# Observed test statistic: difference in means
observed_stat <- mean(group1) - mean(group2)

# Combine all data
combined <- c(group1, group2)
group_labels <- c(rep(1, length(group1)), rep(2, length(group2)))

# Perform permutation test
n_perm <- 10000
perm_stats <- numeric(n_perm)
set.seed(123)  # for reproducibility

for (i in 1:n_perm) {
  permuted <- sample(group_labels)
  perm_group1 <- combined[permuted == 1]
  perm_group2 <- combined[permuted == 2]
  perm_stats[i] <- mean(perm_group1) - mean(perm_group2)
}

# Two-sided p-value
p_value <- mean(abs(perm_stats) >= abs(observed_stat))

# Output result
cat("Observed difference in means (CA - BC):", observed_stat, "\n")
#Observed difference in means (CA - BC): 0.7938554 
cat("Permutation test p-value:", p_value, "\n")
#Permutation test p-value: 0.0566

### Lesions by Sex ###
biome.whale %>% count(Sex, sort = TRUE)
#  Sex          n
#  1 Male       10
#. 2 Female      3
#. 3 Unknown     2


# Filter to only Male and Female
biome_sex <- biome.whale %>% 
  filter(Sex %in% c("Male", "Female")) %>%
  drop_na(percentlesions.total)

# Extract data
lesions <- biome_sex$percentlesions.total
sex <- biome_sex$Sex

# Check group sizes
table(sex)
# Female    Male Unknown 
#.   3      10       0 

# Observed difference in means
group1 <- lesions[sex == "Male"]
group2 <- lesions[sex == "Female"]
observed_stat <- mean(group1) - mean(group2)

# Combine data
combined <- c(group1, group2)
group_labels <- c(rep("Male", length(group1)), rep("Female", length(group2)))

# Permutation test
n_perm <- 10000
perm_stats <- numeric(n_perm)
set.seed(123)

for (i in 1:n_perm) {
  permuted <- sample(group_labels)
  perm_group1 <- combined[permuted == "Male"]
  perm_group2 <- combined[permuted == "Female"]
  perm_stats[i] <- mean(perm_group1) - mean(perm_group2)
}

# Two-sided p-value
p_value <- mean(abs(perm_stats) >= abs(observed_stat))

# Output
cat("Observed difference in means (Male - Female):", observed_stat, "\n")
#Observed difference in means (Male - Female): -0.8172684 
cat("Permutation test p-value:", p_value, "\n")
#Permutation test p-value: 0.0739

# Violin + boxplot
biome.whale$Region.order<- factor(biome.whale$Region, levels = c("CA", "BC"))

lesions_region.plot<- ggplot(biome.whale, aes(x = Region.order, y = percentlesions.total, fill = Region)) +
  geom_violin(trim = FALSE, alpha = 0.5) + coord_cartesian(ylim = c(-1.9, 4.5)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  scale_fill_manual(values=c("#21918c","#fde725")) +
  labs(y = "Percentage of VSA with Lesions", x = "Geographic Location") +
  scale_x_discrete(labels = c("BC" = "British Columbia", "CA" = "California")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())  +
  theme(legend.position = "none")
lesions_region.plot

# Violin + boxplot
lesions_sex.plot<- ggplot(biome_sex, aes(x = Sex, y = percentlesions.total, fill = Sex)) +
  geom_violin(trim = FALSE, alpha = 0.5) + coord_cartesian(ylim = c(-1.9, 4.5)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  scale_fill_manual(values=c("#21918c","#fde725")) +
  labs(y = "Percentage of VSA with Lesions", x = "Sex") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())  +
  theme(legend.position = "none")
lesions_sex.plot

##### Plot both ###
lesions_factor<- ggarrange(lesions_region.plot, lesions_sex.plot, nrow = 1,ncol = 2, 
                           labels = c("(a)", "(b)"), font.label=list(color="black",size=11),
                           common.legend = FALSE,  heights = c(1, 0.5, 1))
lesions_factor

ggsave(
  filename = "FigS13_lesions_by_factor.png",
  plot = lesions_factor,
  width = 10,        # in inches
  height = 5,        # in inches
  units = "in",
  dpi = 800
)

###########################################
## Statistical Analysis - Total ESVs
##########################################
#Model Testing
str(biome.whale)

model1<-lm(Total.ESVs~ Region + Age + Sex + Month + percentlesions.total + scar_PC1 + scar_PC2, data=biome.whale)
summary(model1)
plot(model1) 
performance::check_model(model1) 
AIC(model1)
#[1] 122.2857

model2<-lm(Total.ESVs~ Region + Age + Sex + Month + Number.likely.lesion.types + scar_PC1 + scar_PC2, data=biome.whale)
summary(model2)
performance::check_model(model2)
plot(model2) 
AIC(model2)
#[1] 122.0576

model3<-lm(Total.ESVs~ Region + Age + Sex + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model3)
plot(model3) 
performance::check_model(model3) 
AIC(model3)
#[1] 130.177

model4<-lm(Total.ESVs~ Region + Age + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model4)
plot(model4) 
performance::check_model(model4) 
AIC(model4)
#[1] 129.2374

model5<-lm(Total.ESVs~ Region + Age + Month + percentlesions.total, data=biome.whale)
summary(model5)
plot(model5) 
performance::check_model(model5) 
AIC(model5)
#[1] 128.2051

biome.whale$log_ESVs <- log(biome.whale$Total.ESVs + 1)
model6 <- lm(log_ESVs ~  Region + Age + Month + percentlesions.total, data = biome.whale)
summary(model6)
plot(model6)
performance::check_model(model6)
AIC(model6)
#[1] 28.64816 

model7 <- lm(log_ESVs ~ Age + Month + percentlesions.total, data = biome.whale)
summary(model7)
plot(model7)
performance::check_model(model7)
AIC(model7)
#[1] 12.55393 *best model - final model - explains 66% of variation

#lm(formula = log_ESVs ~ Age + Month + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.5850 -0.2589 -0.0701  0.3569  0.5850 
#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           3.000364   0.320267   9.368 1.38e-05 ***
#Age                  -0.017788   0.006305  -2.821  0.02245 *  
#MonthFebruary         2.132030   0.268141   7.951 4.56e-05 ***
#MonthJanuary          1.409507   0.302941   4.653  0.00164 ** 
#MonthNovember         1.311002   0.280340   4.676  0.00159 ** 
#MonthOctober         -0.313009   0.383682  -0.816  0.43823    
#percentlesions.total -0.121001   0.140135  -0.863  0.41302    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2954 on 8 degrees of freedom
#Multiple R-squared:  0.9524,	Adjusted R-squared:  0.9167 
#F-statistic: 26.67 on 6 and 8 DF,  p-value: 7.13e-05

model8 <- lm(log_ESVs ~ Month + percentlesions.total, data = biome.whale)
summary(model8)
plot(model8)
performance::check_model(model8)
AIC(model8)
#[1] 20.91292

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "January")
model7a <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7a)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.41667 -0.16075  0.04539  0.11110  0.41667 

#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           4.409872   0.200219  22.025 1.91e-08 ***
#Age                  -0.017788   0.006305  -2.821  0.02245 *  
#Month.levelApril     -1.409507   0.302941  -4.653  0.00164 ** 
#Month.levelFebruary   0.722523   0.218738   3.303  0.01081 *  
#Month.levelNovember  -0.098505   0.257969  -0.382  0.71252    
#Month.levelOctober   -1.722516   0.349118  -4.934  0.00114 ** 
#percentlesions.total -0.121001   0.140135  -0.863  0.41302    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2954 on 8 degrees of freedom
#Multiple R-squared:  0.9524,	Adjusted R-squared:  0.9167 
#F-statistic: 26.67 on 6 and 8 DF,  p-value: 7.13e-05

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "February")
model7b <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7b)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.41667 -0.16075  0.04539  0.11110  0.41667 

#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           5.132395   0.191883  26.748 4.11e-09 ***
#Age                  -0.017788   0.006305  -2.821  0.02245 *  
#Month.levelApril     -2.132030   0.268141  -7.951 4.56e-05 ***
#Month.levelJanuary   -0.722523   0.218738  -3.303  0.01081 *  
#Month.levelNovember  -0.821028   0.227733  -3.605  0.00693 ** 
#Month.levelOctober   -2.445039   0.338013  -7.234 8.95e-05 ***
#percentlesions.total -0.121001   0.140135  -0.863  0.41302    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2954 on 8 degrees of freedom
#Multiple R-squared:  0.9524,	Adjusted R-squared:  0.9167 
#F-statistic: 26.67 on 6 and 8 DF,  p-value: 7.13e-05

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "October")
model7c <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7c)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.41667 -0.16075  0.04539  0.11110  0.41667 

#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           2.687356   0.337383   7.965 4.51e-05 ***
#Age                  -0.017788   0.006305  -2.821  0.02245 *  
#Month.levelApril      0.313009   0.383682   0.816  0.43823    
#Month.levelFebruary   2.445039   0.338013   7.234 8.95e-05 ***
#Month.levelJanuary    1.722516   0.349118   4.934  0.00114 ** 
#Month.levelNovember   1.624010   0.323384   5.022  0.00102 ** 
#percentlesions.total -0.121001   0.140135  -0.863  0.41302    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2954 on 8 degrees of freedom
#Multiple R-squared:  0.9524,	Adjusted R-squared:  0.9167 
#F-statistic: 26.67 on 6 and 8 DF,  p-value: 7.13e-05

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "November")
model7d <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7d)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.41667 -0.16075  0.04539  0.11110  0.41667 

#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           4.311366   0.256157  16.831 1.57e-07 ***
#Age                  -0.017788   0.006305  -2.821  0.02245 *  
#Month.levelApril     -1.311002   0.280340  -4.676  0.00159 ** 
#Month.levelFebruary   0.821028   0.227733   3.605  0.00693 ** 
#Month.levelJanuary    0.098505   0.257969   0.382  0.71252    
#Month.levelOctober   -1.624010   0.323384  -5.022  0.00102 ** 
#percentlesions.total -0.121001   0.140135  -0.863  0.41302    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2954 on 8 degrees of freedom
#Multiple R-squared:  0.9524,	Adjusted R-squared:  0.9167 
#F-statistic: 26.67 on 6 and 8 DF,  p-value: 7.13e-05

### Visual best model ###
plot(allEffects(model7a)) #January relevel

# Predict effects for each variable

#Effect of Age on log(ESVs)"
effect_age <- ggpredict(model7a, terms = "Age")

effect_age.plot<- ggplot(effect_age, aes(x = x, y = predicted)) +
  geom_line(color = "steelblue", linewidth = 1) + coord_cartesian(ylim = c(2, 6)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  labs(x = "Age", y = "Predicted log(Total ESVs)") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
effect_age.plot

#Effect of Month on log(ESVs)
effect_month <- ggpredict(model7a, terms = "Month.level")
month_order <- c("January", "February", "April","October", "November") # Define month order
effect_month$x <- factor(effect_month$x, levels = month_order) # Reorder the Month.level variable (x) in model


effect_month.plot<- ggplot(effect_month, aes(x = x, y = predicted, group = 1)) +
  geom_point(size = 3, color = "steelblue") + coord_cartesian(ylim = c(2, 6)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_line(color = "steelblue", linetype = "dashed", size = 1) +  # Trend line
  labs(x = "Month", y = "Predicted log(Total ESVs)") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
effect_month.plot

#Effect of Lesion % on log(ESVs)
effect_lesions <- ggpredict(model7a, terms = "percentlesions.total")

effect_lesions.plot<- ggplot(effect_lesions, aes(x = x, y = predicted)) +
  geom_line(color = "steelblue", size = 1) + coord_cartesian(ylim = c(2, 6)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  labs(x = "Percentage of VSA with Lesions", y = "Predicted log(Total ESVs)") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) 
effect_lesions.plot
##### Plot all effects ###
effects_plot<- ggarrange(effect_age.plot, effect_lesions.plot, effect_month.plot, nrow = 1,ncol = 3, 
                         labels = c("(a)", "(b)", "(c)"), font.label=list(color="black",size=11),
                         common.legend = FALSE,  heights = c(1, 0.5, 1))
effects_plot
ggsave(
  filename = "FigS14_lesions_effect_plot.png",
  plot = effects_plot,
  width = 10,        # in inches
  height = 5,        # in inches
  units = "in",
  dpi = 800
)

#################################################
## Statistical Analysis - Total Bacteria Genera
################################################

model.1<-lm(Total.bacteria.genera~ Region + Age + Sex + Month + percentlesions.total + scar_PC1 + scar_PC2, data=biome.whale)
summary(model.1)
plot(model.1) 
performance::check_model(model.1) 
AIC(model.1)
#[1] 125.9515


model.2<-lm(Total.bacteria.genera~ Region + Age + Sex + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model.2)
plot(model.2) 
performance::check_model(model.2) 
AIC(model.2)
#[1] 124.609


biome.whale$log_Total.bacteria.genera <- log(biome.whale$Total.bacteria.genera + 1)
model.3<-lm(Total.bacteria.genera ~ Region + Sex + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model.3)
plot(model.3) 
performance::check_model(model.3) 
AIC(model.3)
#[1] 122.721


model.4<-lm(Total.bacteria.genera ~ Region + Month + percentlesions.total, data=biome.whale)
summary(model.4)
plot(model.4) 
performance::check_model(model.4) 
AIC(model.4)
#[1] 127.2506


model.5<-lm(Total.bacteria.genera ~ Region + Sex + Month, data=biome.whale)
summary(model.5)
plot(model.5) 
performance::check_model(model.5) 
AIC(model.5)
#[1] 120.6175 - best model, explains 96% of variance


#Residuals:
#Min      1Q  Median      3Q     Max 
#-23.409  -2.000   0.000   3.568  14.591 

#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    120.000     10.833  11.077 3.94e-06 ***
#RegionCA      -110.500     10.833 -10.200 7.32e-06 ***
#SexMale         30.000     13.268   2.261  0.05363 .  
#SexUnknown      19.955     16.002   1.247  0.24766    
#MonthFebruary   35.409      9.239   3.833  0.00500 ** 
#MonthJanuary   -39.318     10.329  -3.807  0.00519 ** 
#MonthNovember       NA         NA      NA       NA    
#MonthOctober    40.000     17.129   2.335  0.04777 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.83 on 8 degrees of freedom
#Multiple R-squared:  0.9797,	Adjusted R-squared:  0.9645 
#F-statistic: 64.39 on 6 and 8 DF,  p-value: 2.46e-06


biome.whale$Month.level <- relevel(biome.whale$Month, ref = "January")
model5a <- lm(Total.bacteria.genera~ Region + Sex + Month.level, data=biome.whale)
summary(model5a)

#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           80.682     14.968   5.390 0.000654 ***
#RegionCA             -31.182     16.814  -1.854 0.100784    
#SexMale               30.000     13.268   2.261 0.053630 .  
#SexUnknown            19.955     16.002   1.247 0.247657    
#Month.levelApril     -40.000     17.129  -2.335 0.047768 *  
#Month.levelFebruary   74.727      8.001   9.340 1.41e-05 ***
#Month.levelNovember   39.318     10.329   3.807 0.005188 ** 
#Month.levelOctober        NA         NA      NA       NA    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.83 on 8 degrees of freedom
#Multiple R-squared:  0.9797,	Adjusted R-squared:  0.9645 
#F-statistic: 64.39 on 6 and 8 DF,  p-value: 2.46e-06


biome.whale$Month.level <- relevel(biome.whale$Month, ref = "February")
model5b <- lm(Total.bacteria.genera~ Region + Sex + Month.level, data=biome.whale)
summary(model5b)

                      #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          155.409     14.238  10.915 4.40e-06 ***
#RegionCA            -105.909     16.168  -6.551 0.000178 ***
#SexMale               30.000     13.268   2.261 0.053630 .  
#SexUnknown            19.955     16.002   1.247 0.247657    
#Month.levelApril     -40.000     17.129  -2.335 0.047768 *  
#Month.levelJanuary   -74.727      8.001  -9.340 1.41e-05 ***
#Month.levelNovember  -35.409      9.239  -3.833 0.004998 ** 
#Month.levelOctober        NA         NA      NA       NA    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.83 on 8 degrees of freedom
#Multiple R-squared:  0.9797,	Adjusted R-squared:  0.9645 
#F-statistic: 64.39 on 6 and 8 DF,  p-value: 2.46e-06


biome.whale$Month.level <- relevel(biome.whale$Month, ref = "April")
model5c <- lm(Total.bacteria.genera~ Region + Sex + Month.level, data=biome.whale)
summary(model5c)

#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          120.000     10.833  11.077 3.94e-06 ***
#RegionCA            -110.500     10.833 -10.200 7.32e-06 ***
#SexMale               30.000     13.268   2.261  0.05363 .  
#SexUnknown            19.955     16.002   1.247  0.24766    
#Month.levelFebruary   35.409      9.239   3.833  0.00500 ** 
#Month.levelJanuary   -39.318     10.329  -3.807  0.00519 ** 
#Month.levelNovember       NA         NA      NA       NA    
#Month.levelOctober    40.000     17.129   2.335  0.04777 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.83 on 8 degrees of freedom
#Multiple R-squared:  0.9797,	Adjusted R-squared:  0.9645 
#F-statistic: 64.39 on 6 and 8 DF,  p-value: 2.46e-06


biome.whale$Month.level <- relevel(biome.whale$Month, ref = "October")
model5d <- lm(Total.bacteria.genera~ Region + Sex + Month.level, data=biome.whale)
summary(model5d)

#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          120.000     10.833  11.077 3.94e-06 ***
#RegionCA             -70.500     13.268  -5.314 0.000717 ***
#SexMale               30.000     13.268   2.261 0.053630 .  
#SexUnknown            19.955     16.002   1.247 0.247657    
#Month.levelApril     -40.000     17.129  -2.335 0.047768 *  
#Month.levelFebruary   35.409      9.239   3.833 0.004998 ** 
#Month.levelJanuary   -39.318     10.329  -3.807 0.005188 ** 
#Month.levelNovember       NA         NA      NA       NA    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.83 on 8 degrees of freedom
#Multiple R-squared:  0.9797,	Adjusted R-squared:  0.9645 
#F-statistic: 64.39 on 6 and 8 DF,  p-value: 2.46e-06




biome.whale$Month.level <- relevel(biome.whale$Month, ref = "November")
model5e <- lm(Total.bacteria.genera~ Region + Sex + Month.level, data=biome.whale)
summary(model5e)

                  #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          120.000     10.833  11.077 3.94e-06 ***
#RegionCA             -70.500     13.268  -5.314 0.000717 ***
#SexMale               30.000     13.268   2.261 0.053630 .  
#SexUnknown            19.955     16.002   1.247 0.247657    
#Month.levelApril     -40.000     17.129  -2.335 0.047768 *  
#Month.levelFebruary   35.409      9.239   3.833 0.004998 ** 
#Month.levelJanuary   -39.318     10.329  -3.807 0.005188 ** 
#Month.levelOctober        NA         NA      NA       NA    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.83 on 8 degrees of freedom
#Multiple R-squared:  0.9797,	Adjusted R-squared:  0.9645 
#F-statistic: 64.39 on 6 and 8 DF,  p-value: 2.46e-06




### Visual best model ###
plot(allEffects(model5a)) #January relevel

#Effect of sex on log(genera)
effect_sex2 <- ggpredict(model5a, terms = "Sex") 

effect_sex.plot2 <- ggplot(effect_sex2, aes(x = x, y = predicted, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  coord_cartesian(ylim = c(-50, 130)) +
  labs(x = "Sex", y = "Predicted log(Total Bacteria Genera)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

effect_sex.plot2

#Effect of Month on log(genera)
effect_month2 <- ggpredict(model5a, terms = "Month.level") 
month_order <- c("January", "February", "April","October", "November") # Define month order
effect_month2$x <- factor(effect_month2$x, levels = month_order) # Reorder the Month.level variable (x) in model

effect_month.plot2<- ggplot(effect_month2, aes(x = x, y = predicted, group = 1)) +
  geom_point(size = 3, color = "steelblue") + coord_cartesian(ylim = c(-50, 190)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_line(color = "steelblue", linetype = "dashed", linewidth = 1) +  # Trend line
  labs(x = "Month", y = "Predicted log(Total Bacteria Genera)") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
effect_month.plot2

#Effect of region on log(genera)
effect_region2 <- ggpredict(model5a, terms = "Region")

effect_region.plot2 <- ggplot(effect_region2, aes(x = x, y = predicted, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  coord_cartesian(ylim = c(-50, 130)) +
  labs(x = "Region", y = "Predicted log(Total Bacteria Genera)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

effect_region.plot2

##### Plot all effects ###
effects_plot2<- ggarrange(effect_region.plot2, effect_sex.plot2, effect_month.plot2, nrow = 1,ncol = 3, 
                          labels = c("(a)", "(b)", "(c)"), font.label=list(color="black",size=11),
                          common.legend = FALSE,  heights = c(1, 0.5, 1))
effects_plot2

ggsave(
  filename = "FigS15_genera_effect_plot.png",
  plot = effects_plot2,
  width = 10,        # in inches
  height = 5,        # in inches
  units = "in",
  dpi = 800
)
##### END #####