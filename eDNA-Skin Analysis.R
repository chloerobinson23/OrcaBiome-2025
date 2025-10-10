### Orcabiome eDNA-Skin Models ###
#Created by: Brittany Visona-Kelly
#Created: June 6, 2025

###########################
## Data Processing
###########################

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

#install.packages("DescTools")
library(DescTools)
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
#library(dplyr)
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
#install.packages("GGally")
#library(GGally)
GGally::ggpairs(biome.whale[, c("Age", "Number.likely.lesion.types","percentlesions.total", "scar_PC1", "scar_PC2")])
#Number.likely.lesion.types and percentlesions.total are correlated (> 0.80). Only use 1 in models. 

#Test normality
#install.packages("ggpubr")
#library(ggpubr)
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

library(tidyr)
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
#ggexport(lesions_factor, filename = "FigS9_lesions by factor.pdf", width = 10, height = 5, units = "in", dpi = 800)

###########################################
## Statistical Analysis - Total ESVs
##########################################
#Model Testing
#install.packages("performance")
library(performance)
str(biome.whale)

model1<-lm(Total.ESVs~ Region + Age + Sex + Month + percentlesions.total + scar_PC1 + scar_PC2, data=biome.whale)
summary(model1)
plot(model1) 
performance::check_model(model1) 
AIC(model1)
#[1] 175.0729

model2<-lm(Total.ESVs~ Region + Age + Sex + Month + Number.likely.lesion.types + scar_PC1 + scar_PC2, data=biome.whale)
summary(model2)
performance::check_model(model2)
plot(model2) 
AIC(model2)
#[1] 178.6724

model3<-lm(Total.ESVs~ Region + Age + Sex + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model3)
plot(model3) 
performance::check_model(model3) 
AIC(model3)
#[1] 174.251

model4<-lm(Total.ESVs~ Region + Age + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model4)
plot(model4) 
performance::check_model(model4) 
AIC(model4)
#[1] 171.3593

model5<-lm(Total.ESVs~ Region + Age + Month + percentlesions.total, data=biome.whale)
summary(model5)
plot(model5) 
performance::check_model(model5) 
AIC(model5)
#[1] 169.9911

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
#[1] 28.64816 *best model - final model - explains 66% of variation

#lm(formula = log_ESVs ~ Age + Month + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.5850 -0.2589 -0.0701  0.3569  0.5850 
#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           5.35885    0.54765   9.785 9.98e-06 ***
#  Age                  -0.01682    0.01078  -1.560   0.1574    
#MonthFebruary        -0.53102    0.45851  -1.158   0.2802    
#MonthJanuary         -1.12047    0.51802  -2.163   0.0625 .  
#MonthNovember         0.10554    0.47937   0.220   0.8313    
#MonthOctober         -1.71004    0.65608  -2.606   0.0313 *  
#  percentlesions.total  0.16507    0.23963   0.689   0.5104    
#Residual standard error: 0.5051 on 8 degrees of freedom
#Multiple R-squared:  0.662,	Adjusted R-squared:  0.4085 
#F-statistic: 2.611 on 6 and 8 DF,  p-value: 0.1048

model8 <- lm(log_ESVs ~ Month + percentlesions.total, data = biome.whale)
summary(model8)
plot(model8)
performance::check_model(model8)
AIC(model8)
#[1] 30.63118

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "January")
model7a <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7a)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.5850 -0.2589 -0.0701  0.3569  0.5850 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           4.23839    0.34237  12.380 1.69e-06 ***
#  Age                  -0.01682    0.01078  -1.560   0.1574    
#Month.levelApril      1.12047    0.51802   2.163   0.0625   
#Month.levelFebruary   0.58945    0.37404   1.576   0.1537    
#Month.levelNovember   1.22600    0.44112   2.779   0.0240 *  
#  Month.levelOctober   -0.58957    0.59698  -0.988   0.3523    
#percentlesions.total  0.16507    0.23963   0.689   0.5104    
#Residual standard error: 0.5051 on 8 degrees of freedom
#Multiple R-squared:  0.662,	Adjusted R-squared:  0.4085 
#F-statistic: 2.611 on 6 and 8 DF,  p-value: 0.1048

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "February")
model7b <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7b)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.5850 -0.2589 -0.0701  0.3569  0.5850 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           4.82784    0.32811  14.714 4.47e-07 ***
#  Age                  -0.01682    0.01078  -1.560   0.1574    
#Month.levelApril      0.53102    0.45851   1.158   0.2802    
#Month.levelJanuary   -0.58945    0.37404  -1.576   0.1537    
#Month.levelNovember   0.63655    0.38942   1.635   0.1408    
#Month.levelOctober   -1.17902    0.57799  -2.040   0.0757 
#percentlesions.total  0.16507    0.23963   0.689   0.5104    
#Residual standard error: 0.5051 on 8 degrees of freedom
#Multiple R-squared:  0.662,	Adjusted R-squared:  0.4085 
#F-statistic: 2.611 on 6 and 8 DF,  p-value: 0.1048

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "October")
model7c <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7c)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.5850 -0.2589 -0.0701  0.3569  0.5850 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           3.64881    0.57691   6.325 0.000227 ***
#  Age                  -0.01682    0.01078  -1.560 0.157423    
#Month.levelApril      1.71004    0.65608   2.606 0.031303 *  
#  Month.levelFebruary   1.17902    0.57799   2.040 0.075690 .  
#Month.levelJanuary    0.58957    0.59698   0.988 0.352277    
#Month.levelNovember   1.81558    0.55298   3.283 0.011133 *  
#  percentlesions.total  0.16507    0.23963   0.689 0.510389    
#Residual standard error: 0.5051 on 8 degrees of freedom
#Multiple R-squared:  0.662,	Adjusted R-squared:  0.4085 
#F-statistic: 2.611 on 6 and 8 DF,  p-value: 0.1048

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "November")
model7d <- lm(log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model7d)
#lm(formula = log_ESVs ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.5850 -0.2589 -0.0701  0.3569  0.5850 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           5.46439    0.43802  12.475 1.59e-06 ***
#  Age                  -0.01682    0.01078  -1.560   0.1574    
#Month.levelApril     -0.10554    0.47937  -0.220   0.8313    
#Month.levelFebruary  -0.63655    0.38942  -1.635   0.1408    
#Month.levelJanuary   -1.22600    0.44112  -2.779   0.0240 *  
#  Month.levelOctober   -1.81558    0.55298  -3.283   0.0111 *  
#  percentlesions.total  0.16507    0.23963   0.689   0.5104    
#Residual standard error: 0.5051 on 8 degrees of freedom
#Multiple R-squared:  0.662,	Adjusted R-squared:  0.4085 
#F-statistic: 2.611 on 6 and 8 DF,  p-value: 0.1048

### Visual best model ###
#install.packages("effects")  
library(effects)
plot(allEffects(model7a)) #January relevel

# Predict effects for each variable
#install.packages("ggeffects")
library(ggeffects)

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
#ggexport(effects_plot, filename = "FigS10_effects plot ESV model.pdf", width = 10, height = 5, units = "in", dpi = 800)

#################################################
## Statistical Analysis - Total Bacteria Genera
################################################

model.1<-lm(Total.bacteria.genera~ Region + Age + Sex + Month + percentlesions.total + scar_PC1 + scar_PC2, data=biome.whale)
summary(model.1)
plot(model.1) 
performance::check_model(model.1) 
AIC(model.1)
#[1] 146.2591


model.2<-lm(Total.bacteria.genera~ Region + Age + Sex + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model.2)
plot(model.2) 
performance::check_model(model.2) 
AIC(model.2)
#[1] 144.402

biome.whale$log_Total.bacteria.genera <- log(biome.whale$Total.bacteria.genera + 1)
model.3<-lm(Total.bacteria.genera ~ Age + Month + percentlesions.total + scar_PC1, data=biome.whale)
summary(model.3)
plot(model.3) 
performance::check_model(model.3) 
AIC(model.3)
#[1] 141.6886

model.4<-lm(Total.bacteria.genera ~ Age + Month + percentlesions.total, data=biome.whale)
summary(model.4)
plot(model.4) 
performance::check_model(model.4) 
AIC(model.4)
#[1] 139.9079 *best model - final model - explains 72% of variation
#lm(formula = Total.bacteria.genera ~ Age + Month + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-23.657 -12.822   1.343  13.355  24.615 
#Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)  
#(Intercept)           60.4127    22.3435   2.704   0.0269 *
#  Age                   -0.6774     0.4399  -1.540   0.1621  
#MonthFebruary         17.5917    18.7069   0.940   0.3745  
#MonthJanuary         -10.9299    21.1347  -0.517   0.6190  
#MonthNovember         35.9662    19.5580   1.839   0.1032  
#MonthOctober         -54.3817    26.7676  -2.032   0.0767 .
#percentlesions.total  14.2292     9.7765   1.455   0.1836  
#Residual standard error: 20.61 on 8 degrees of freedom
#Multiple R-squared:  0.7204,	Adjusted R-squared:  0.5108 
#F-statistic: 3.436 on 6 and 8 DF,  p-value: 0.05542

model.5<-lm(Total.bacteria.genera ~ Month + percentlesions.total, data=biome.whale)
summary(model.5)
plot(model.5) 
performance::check_model(model.5) 
AIC(model.5)
#[1] 141.8024

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "January")
model4a <- lm(Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model4a)
#lm(formula = Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-23.657 -12.822   1.343  13.355  24.615 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)           49.4827    13.9683   3.543  0.00759 **
#  Age                   -0.6774     0.4399  -1.540  0.16213   
#Month.levelApril      10.9299    21.1347   0.517  0.61904   
#Month.levelFebruary   28.5216    15.2603   1.869  0.09856 . 
#Month.levelNovember   46.8961    17.9972   2.606  0.03134 * 
#  Month.levelOctober   -43.4518    24.3563  -1.784  0.11226   
#percentlesions.total  14.2292     9.7765   1.455  0.18364   
#Residual standard error: 20.61 on 8 degrees of freedom
#Multiple R-squared:  0.7204,	Adjusted R-squared:  0.5108 
#F-statistic: 3.436 on 6 and 8 DF,  p-value: 0.05542

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "February")
model4b <- lm(Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model4b)
#lm(formula = Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-23.657 -12.822   1.343  13.355  24.615 
#oefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           78.0044    13.3867   5.827 0.000393 ***
#  Age                   -0.6774     0.4399  -1.540 0.162126    
#Month.levelApril     -17.5917    18.7069  -0.940 0.374543    
#Month.levelJanuary   -28.5216    15.2603  -1.869 0.098557 .  
#Month.levelNovember   18.3745    15.8878   1.157 0.280836    
#Month.levelOctober   -71.9734    23.5815  -3.052 0.015770 *  
#  percentlesions.total  14.2292     9.7765   1.455 0.183636    
#Residual standard error: 20.61 on 8 degrees of freedom
#Multiple R-squared:  0.7204,	Adjusted R-squared:  0.5108 
#F-statistic: 3.436 on 6 and 8 DF,  p-value: 0.05542

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "October")
model4c <- lm(Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model4c)
#lm(formula = Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-23.657 -12.822   1.343  13.355  24.615 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)            6.0310    23.5375   0.256  0.80424   
#Age                   -0.6774     0.4399  -1.540  0.16213   
#Month.levelApril      54.3817    26.7676   2.032  0.07666 . 
#Month.levelFebruary   71.9734    23.5815   3.052  0.01577 * 
#  Month.levelJanuary    43.4518    24.3563   1.784  0.11226   
#Month.levelNovember   90.3479    22.5609   4.005  0.00392 **
#  percentlesions.total  14.2292     9.7765   1.455  0.18364   
#Residual standard error: 20.61 on 8 degrees of freedom
#Multiple R-squared:  0.7204,	Adjusted R-squared:  0.5108 
#F-statistic: 3.436 on 6 and 8 DF,  p-value: 0.05542

biome.whale$Month.level <- relevel(biome.whale$Month, ref = "November")
model4d <- lm(Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
summary(model4d)
#lm(formula = Total.bacteria.genera ~ Age + Month.level + percentlesions.total, data = biome.whale)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-23.657 -12.822   1.343  13.355  24.615 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           96.3788    17.8708   5.393 0.000651 ***
#  Age                   -0.6774     0.4399  -1.540 0.162126    
#Month.levelApril     -35.9662    19.5580  -1.839 0.103211    
#Month.levelFebruary  -18.3745    15.8878  -1.157 0.280836    
#Month.levelJanuary   -46.8961    17.9972  -2.606 0.031337 *  
#  Month.levelOctober   -90.3479    22.5609  -4.005 0.003924 ** 
#  percentlesions.total  14.2292     9.7765   1.455 0.183636    
#Residual standard error: 20.61 on 8 degrees of freedom
#Multiple R-squared:  0.7204,	Adjusted R-squared:  0.5108 
#F-statistic: 3.436 on 6 and 8 DF,  p-value: 0.05542

### Visual best model ###
plot(allEffects(model4a)) #January relevel

#Effect of Age on log(ESVs)"
effect_age2 <- ggpredict(model4a, terms = "Age")

effect_age.plot2<- ggplot(effect_age2, aes(x = x, y = predicted)) +
  geom_line(color = "steelblue", linewidth = 1) + coord_cartesian(ylim = c(-50, 130)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  labs(x = "Age", y = "Predicted log(Total Bacteria Genera)") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
effect_age.plot2

#Effect of Month on log(ESVs)
effect_month2 <- ggpredict(model4a, terms = "Month.level") 
month_order <- c("January", "February", "April","October", "November") # Define month order
effect_month2$x <- factor(effect_month2$x, levels = month_order) # Reorder the Month.level variable (x) in model

effect_month.plot2<- ggplot(effect_month2, aes(x = x, y = predicted, group = 1)) +
  geom_point(size = 3, color = "steelblue") + coord_cartesian(ylim = c(-50, 130)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_line(color = "steelblue", linetype = "dashed", size = 1) +  # Trend line
  labs(x = "Month", y = "Predicted log(Total Bacteria Genera)") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
effect_month.plot2

#Effect of Lesion % on log(ESVs)
effect_lesions2 <- ggpredict(model4a, terms = "percentlesions.total")

effect_lesions.plot2<- ggplot(effect_lesions2, aes(x = x, y = predicted)) +
  geom_line(color = "steelblue", size = 1) + coord_cartesian(ylim = c(-50, 130)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  labs(x = "Percentage of VSA with Lesions", y = "Predicted log(Total Bacteria Genera)") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) 
effect_lesions.plot2

##### Plot all effects ###
effects_plot2<- ggarrange(effect_age.plot2, effect_lesions.plot2, effect_month.plot2, nrow = 1,ncol = 3, 
                          labels = c("(a)", "(b)", "(c)"), font.label=list(color="black",size=11),
                          common.legend = FALSE,  heights = c(1, 0.5, 1))
effects_plot2
#ggexport(effects_plot2, filename = "FigS11_effects plot bacteria genera model.pdf", width = 10, height = 5, units = "in", dpi = 800)


### END ###