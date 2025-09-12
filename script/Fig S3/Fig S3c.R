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

setwd("G:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/")

# 1. Taxonomy Table
taxa <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/tax_final250_220.rds")
# check your taxonomic classifications #
taxa.print<- taxa
rownames(taxa.print)
head(taxa.print)
unname(head(taxa))

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

###ASV table G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/
asv.table<- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/seqtab_final250_220.rds")
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)


#Now we can make the phyloseq object
ps2 <- phyloseq(asv.table2, tax_table(taxa), sample_data(meta))

# while it doesn't seem to be the case, if it were we would use the following lines to remove samples with low sequence counts
set.seed(500)
ps.2 = subset_taxa(ps2, Kingdom == "Bacteria" | Kingdom == "Archaea")
ps.pruned <- prune_samples(sample_sums(ps.2)>=2, ps.2)
ps.3 = subset_samples(ps.pruned, Compartment !=  "BLANK")
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))

ps.perc.3 <- transform_sample_counts(ps.3.SC, function(x) x / sum(x)) 
ps.pruned_taxa <- filter_taxa(ps.perc.3, function(x) sum(x) > .005, TRUE)

#####PCoA for ASV-level data with Bray-Curtis
PCoA <- ordinate(ps.perc.3, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x,ps.perc.3.H@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
nrow(x)
'bronw'
cols <- c("Frond" = "#117733", "Root" = "#8B4513", "Water" = "#87CEFA", "Inoculum" = "black")

x$Compartment[which(x$Compartment == 'Front')] <- "Frond"
x2 = x[x$Compartment!='Water',]
unique(x$Specie.1)



ggplot(x, aes(Axis.1, Axis.2)) + #
  geom_point(aes(color=Compartment), size=3) + #
  theme_few() +
  theme(legend.position="top")+ 
  scale_colour_manual(values = cols) +
  labs(y=ylab_text, x=xlab_text)

