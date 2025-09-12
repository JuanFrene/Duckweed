packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')

sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

# Control sequence ####
track <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/track.rds")
write.table(track, "track 16S.txt")

# 1. Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

###ASV table
asv.table <- readRDS('G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/seqtab_final.rds')
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

##ASing to the syncom
taxaSyncom <- assignTaxonomy(asv.table2, "G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/Syncom_Sequence.pcr.unique.fasta", multithread=TRUE)#, minBoot = 98

#write.table(taxaSyncom, "taxaSyncom.txt")
#taxaSyncom = read.table("taxaSyncom.txt", header = TRUE, row.names = 1)

#Now we can make the phyloseq object
ps.Syncom <- phyloseq(asv.table2, tax_table(taxaSyncom), sample_data(meta))

taxa_names(ps.Syncom) <- paste0("ASV", seq(ntaxa(ps.Syncom)))

ps.Syncom = subset_taxa(ps.Syncom, Genus  != "NA")
#ps.SyncomSilva = subset_taxa(ps.SyncomSilva, Genus  != "NA")

ps.Syncom@sam_data
Tax95s=ps.Syncom@tax_table
Tax95s=data.frame(Tax95s)
unique(Tax95s$Kingdom)

ps.Syncom.3 = subset_samples(ps.Syncom, ID !=  "C-")
ps.pruned <- prune_taxa(taxa_sums(ps.Syncom.3)>=1, ps.Syncom.3)
ps.pruned3 <- prune_samples(sample_sums(ps.pruned)>1, ps.pruned) 
ps.perc <- transform_sample_counts(ps.pruned3, function(x) x / sum(x)) 


ps.Syncom2 = subset_samples(ps.Syncom, Species !=  "LP8539")
ps.Syncom3 = subset_samples(ps.Syncom2, Species !=  "LY9250")

ps.perc.Control = subset_samples(ps.Syncom3, Drugs ==  "Control")
ps.perc.DMSO = subset_samples(ps.Syncom, Drugs ==  "DMSO")



ps.perc.Lysine = subset_samples(ps.Syncom3, Drugs ==  "Lysine")
ps.perc.Lysine2 <- prune_samples(sample_sums(ps.perc.Lysine)>1, ps.perc.Lysine) 

#####PCoA for ASV-level data with Bray-Curtis
# Calculates Bray-Curtis distances between samples. Because taxa is in
# columns, it is used to compare different samples. We transpose the
# assay to get taxa to columns
PCoA <- ordinate(ps.perc.Lysine2, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x, ps.perc@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Species*Drugs,x,mean)
PCoAmeands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Species*Drugs,x,se)
PCoAmean2 =merge(PCoAmean,PCoAmeands[,-(2)],by="Species")

gg2 <- merge(PCoAmean2,x[,c(9,12)], by="Species")
gg2$Complexity = as.numeric(gg2$Complexity)

ggplot(data = gg2, aes(mean.x, mean.y, label=Species)) +
    geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
    geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
    geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.1, alpha = 0.15) +
    geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),size = 0.1, alpha = 0.15)  + 
    geom_point( size = 6, aes(col = Complexity),stroke = 1) + 
    scale_colour_gradientn(colours=c("Blue","Red"))+
    labs(y=ylab_text, x=xlab_text)+ theme_few()
    
###PERMANOVA
# Calculate bray curtis distance matrix
ps.3_bray <- phyloseq::distance(ps.perc.Citruline2, method = "bray")
# make a data frame from the sample_data
sample_ps.3c <- data.frame(sample_data(ps.perc.Citruline2))
# Adonis test
adonis2(ps.3_bray ~ Species, data = sample_ps.3c)


#Lysine: 
#adonis2(formula = ps.3_bray ~ Complexity * Species, data = sample_ps.3c)
#Complexity  1   0.5616 0.14202 2.2936  0.060 .
#Species     2   0.6995 0.17688 1.4283  0.144  
#Residual   11   2.6936 0.68110                

