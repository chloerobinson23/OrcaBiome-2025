### OrcaBiome Paper Analysis
### Rarefaction curve
### March 10, 2025


# Load libraries
library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rarecurve
library(purrr) # for map_dfr
library(ggplot2) # ggplot
library(scales) # comma
library(gridExtra) # grid.arrange
library(cowplot) # get_legend

###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#####################################################################

# Read in data
A <- read.csv(file="OrcaBiome_results_0.8_cut_genera_removed_noNTC_v2.csv", head=TRUE)


#######################################################
# Create dataframes for vegan based on all taxa
######################################################

# Split up SampleName with pkg 'stringr'
A.1 <- data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[25:29] <- c("Sample", "SampleType","WhaleID","Sex","Population")
#Sample type = blank or flukeprint, whaleID = ID of whale, sex = sex of whale, population = CA/BC
#edit numbers in [] to correspond with what cols have the above attribute headings

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
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
# 15% = 18618.15

# set random seed
set.seed(1234)

# Do rarefection with pkg 'vegan'
rarecurveout <- rarecurve2(esv.notnull2, 
                           sample=esv.percentile, 
                           step=500, 
                           label=T)

# Reformat vegan list as df (cols OTU, raw.read)
rare.df <- lapply(rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Add sample names to vegan output (df) (rownames)
sample_names <- rownames(esv.notnull2)
names(rare.df) <- sample_names

# Map rownames to vegan output (df)
rare.df <- map_dfr(rare.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
rare.df <- data.frame(rare.df, do.call(rbind, str_split(rare.df$sample,"_")), stringsAsFactors = FALSE)
names(rare.df)[4:8]<-c("Sample","SampleType","WhaleID","Sex","Population")

# Create factors
rare.df$SampleType <- factor(rare.df$SampleType, 
                       levels = c("BL","FP"),
                       labels = c("Seawater Control", "Flukeprint"))
rare.df$Sex <- factor(rare.df$Sex,
                             levels = c("M", "F","U","NA"),
                             labels = c("Male", "Female","Unknown","NA"))
rare.df$Population <- factor(rare.df$Population,
                          levels = c("MB", "BC"),
                          labels = c("California", "British Columbia"))

# Create plot 1 (color by SampleType)
p1.tmp <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Sample Type") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = SampleType), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  scale_color_manual(values = c("Seawater Control" = "lightblue", "Flukeprint" = "black")) +  # Custom colors
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank()) +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 2))) 

l1 <- get_legend(p1.tmp)

p1 <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Sample Type") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = SampleType), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  scale_color_manual(values = c("Seawater Control" = "lightblue", "Flukeprint" = "black")) +  # Custom colors
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_text(),
        legend.position = "none")

p1

# Create plot 2 (color by Population)
p2.tmp <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Population") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Population), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  scale_color_manual(values = c("California" = "#3B4D93", "British Columbia" = "#2A9D8F")) +  # Custom colors
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank()) +
  guides(color=guide_legend(ncol=2, nrow=1, override.aes = list(size = 2)))  # Adjust legend to one row

l2 <- get_legend(p2.tmp)

p2 <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Population") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Population), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  scale_color_manual(values = c("California" = "#3B4D93", "British Columbia" = "#2A9D8F")) +  # Custom colors
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none")

# Arrange the plots and legends in a grid layout
g <- plot_grid(p1, p2, l1, l2, ncol = 2, nrow = 2, rel_heights = c(1, .2))

# Display the final plot
g

ggsave("FigS1_rarefaction.jpeg", plot = g, dpi = 300, height = 6, width = 8, units = "in")
# vertical line represents 15% percentile of read depth across sample (excluding controls)
# rarefaction curve with NTC and PC removed

###### END #######