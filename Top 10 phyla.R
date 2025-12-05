### OrcaBiome Paper Analysis
### Top 10
### December 02, 2025
### Based off Teresita M. Porter, Dec. 23, 2019


# Load libraries
library(reshape2) # dcast
library(ggplot2) # ggplot


# Read in data
A <- read.csv(file="OrcaBiome_results_0.8_cut_NTC ESVs removed_OB36_removed.csv", head=TRUE)

# Summarize ESVs in all unique phyla (pooled across samples)
A.esv <- dcast(A, Phylum ~ . , value.var="GlobalESV", function (x) {length(unique(x))})
names(A.esv)<-c("Phyla","GlobalESV")

# Sort by descending ESVs
A.esv.desc<-A.esv[order(-A.esv$GlobalESV),]

# Summarize reads in all detected phyla
A.read <- dcast(A, Phylum ~ . , value.var="ESVsize", fun.aggregate=sum)
names(A.read)<-c("Phyla","ESVsize")

# Sort by descending reads
A.read.desc<-A.read[order(-A.read$ESVsize),]

# Calc proportions
A.esv.phyla<-A.esv.desc[,1]
A.esvprop<-round(A.esv.desc[,2]/sum(A.esv.desc[,2])*100,digits=2)
A.read.phyla<-A.read.desc[,1]
A.readprop<-round(A.read.desc[,2]/sum(A.read.desc[,2])*100,digits=2)

# Create esv prop table
A.esv.table<-data.frame(phylum=A.esv.phyla, esv=A.esvprop)

# Create read prop table
A.read.table<-data.frame(phylum=A.read.phyla, read=A.readprop)

# Keep top 11 (as have to remove 1 from each due to lack of ESVs/reads)
A.esv.top<-A.esv.table[1:11,]
A.read.top<-A.read.table[1:11,]

# Create 'other' df
A.esv.other<-A.esv.table[-(1:11),]
A.read.other<-A.read.table[-(1:11),]

# Sum 'other' line
A.esv.other.sum <- sum(A.esv.other$esv)
A.read.other.sum <- sum(A.read.other$read)

# Create df record for other
A.esv.other.sum.df <- data.frame("phylum"="Other", "esv"=A.esv.other.sum)
A.read.other.sum.df <- data.frame("phylum"="Other", "read"=A.read.other.sum)

# Add other to top df
A.esv.rbind <- rbind(A.esv.top, A.esv.other.sum.df)
A.read.rbind <- rbind(A.read.top, A.read.other.sum.df)

# Outer join esvs and reads
A.table<-merge(A.esv.rbind, A.read.rbind, by="phylum", all=TRUE)

# Create long form form ggplot
A.long<-melt(A.table, id=c("phylum"))

# Exclude rows which do not have have ESVs vs reads
A.long.2 <- A.long[-c(1, 11, 14, 24),]

# Create factors
A.long.2$variable <- factor(A.long.2$variable,
                          levels=c("esv","read"),
                          labels=c("ESVs","Reads"))

# Convert phylum to factor and reorder levels so "Other" is last
A.long.2$phylum <- factor(A.long.2$phylum, 
                          levels = c(setdiff(unique(A.long.2$phylum), "Other"), "Other"))

# Proprtions for top 10
p<-ggplot(data=A.long.2, aes(x=variable, y=value, fill=phylum, label=value)) +
  geom_bar(stat="identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  labs(y="Proportions", x="Top 10 phyla") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=12),
        axis.title.x=element_blank(),
        plot.title=element_text(size=12),
        legend.position="bottom",
        legend.text=element_text(size=10),
        legend.title=element_blank(),
        legend.key.size = unit(0.45, "cm")) +
  guides(fill = guide_legend(nrow = 3))

p

ggsave("FigS5_Top10Phyla.png", plot = p, dpi = 300, height = 10, width = 8, units = "in")
# NTC and PC controls excluded
# based on raw data filtered at 0.8 genus level

###### END #######

