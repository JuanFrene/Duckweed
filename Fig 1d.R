packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest",'dada2', "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", "reshape2")
sapply(packages, require, character.only = TRUE)              

library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(harrietr)
library(egg)
library(paletteer)
library(ggtree)
library(scales)
library(car)
library(Rmisc)

setwd("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/")

###Track sequencing abundance
#track <- readRDS("C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/track240-200.rds")
#write.table(track, "track 16S 240-240 Exp1.txt")

# 1. Taxonomy Table
taxa <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/tax_final250_220.rds")

# 2. Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

####ASV table
asv.table<- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/seqtab_final250_220.rds")
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

#Now we can make the phyloseq object
ps2 <- phyloseq(asv.table2, tax_table(taxa), sample_data(meta))

ps.2 = subset_taxa(ps2, Kingdom == "Bacteria" | Kingdom == "Archaea")
ps.pruned <- prune_samples(sample_sums(ps.2)>=2, ps.2)
ps.pruned_taxa <- filter_taxa(ps.2, function(x) sum(x) > .005, TRUE)

ps.3 = subset_samples(ps.pruned, Compartment !=  "BLANK")
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))

ps.perc.3 <- transform_sample_counts(ps.3, function(x) x / sum(x)) 


ps.3.Root = subset_samples(ps.3, Compartment ==  "Root")
ps.3.Water = subset_samples(ps.3, Compartment ==  "Water")
ps.3.SInoculum = subset_samples(ps.3, Compartment !=  "Inoculum")
ps.3.SWater = subset_samples(ps.3.SInoculum, Compartment !=  "Water")
ps.3.Front = subset_samples(ps.3, Compartment ==  "Front")
ps.3.SC = subset_samples(ps.3, Specie !=  "Control")
ps.3.SCSI = subset_samples(ps.3.SC, Specie !=  "Inoculum")

ps.3.SWater.perc<- transform_sample_counts(ps.3.SWater, function(x) x / sum(x))
ps.perc.Root = subset_samples(ps.3.SWater.perc, Compartment ==  "Root")
ps.perc.Front = subset_samples(ps.3.SWater.perc, Compartment ==  "Front")

##Transformation
ps.3.SInoculum.H <- transform(ps.3.SInoculum, "hellinger")
ps.3.Root.H <- transform(ps.3.Root, "hellinger")
ps.3.Front.H <- transform(ps.3.Front, "hellinger")

ps.3.Root.H2 = subset_samples(ps.3.Root.H, ID !=  "R58")
ps.3.Root.H3 = subset_samples(ps.3.Root.H2, ID !=  "R65")
ps.3.Root.H4 = subset_samples(ps.3.Root.H3, ID !=  "R65_(2)")


#####PCoA for ASV-level data with Bray-Curtis
# Calculates Bray-Curtis distances between samples. Because taxa is in
# columns, it is used to compare different samples. We transpose the
# assay to get taxa to columns
PCoA <- ordinate(ps.3.Front.H, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x,ps.3.Front.H@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
nrow(x)

PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Specie.1,x,mean)
PCoAds = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Specie.1,x,sd)
PCoAmean <- merge(PCoAds,PCoAmean, by="Specie.1")

gg <- merge(x,PCoAmean,by="Specie.1")

mortality = c(68.85245902,9.677419355,14.94252874,34.375,21.05263158,8.421052632,14.45783133,63.51351351,63.38028169,39.63963964,45,5.376344086,0,1.315789474,1.282051282,12,7.291666667,4.651162791,15.55555556)
RootComp = c(5,5.33,5,6.4,5,5.2,6,5,5.2,6.2,6,7.6,8.4,7.6,6.6,7.25,6.8,6.6,7)
RootComp2 = c(5,5,5.75,6.75,5,5,5.2,5.25,5,5.4,5.6,5.8,7.6,7.333333333,6,6.8,5.8,6,6.6)
Epidermis = c(16.77,14.41,15.3,10.44,18.79,8.07,10.05,14.31,17.55,12.54,7,11.49,32.84,8.06,10.89,10.45,10.44,8.37,9.17)
EpidermisPorc = c(51.86,23.27,17.73,22.56,19.92,23.66,48.81,21.44,26.83,27.71,26.04,18.71,32.88,24.19,41.42,40.69,38.31,39.09,32.25)
Nroot = c(2,1.478333333,1.285714286,1,1,1.946666667,2.266666667,1,1,1,1,1,5.8,2.893333333,3.12,2.406666667,2.82,2.341666667,3.473333333)
PCoAmean = cbind(PCoAmean, mortality,RootComp,RootComp2,Epidermis,EpidermisPorc,Nroot)

ggplot(data = PCoAmean, aes(mean.x, mean.y, )) +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.2, alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),size = 0.2, alpha = 0.4)  + 
  geom_point(shape = 21, size = 6, aes(fill = RootComp),stroke = 1) + 
  #scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42") +
  labs(y=ylab_text, x=xlab_text)+
  scale_fill_gradientn(colours=c('blue','red'))+theme_few() #

