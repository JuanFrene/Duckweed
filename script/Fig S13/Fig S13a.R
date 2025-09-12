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


ps.perc.Methycitrate = subset_samples(ps.Syncom3, Drugs ==  "Methycitrate")
ps.perc.Methycitrate2 <- prune_samples(sample_sums(ps.perc.Methycitrate)>1, ps.perc.Methycitrate) 

ps.perc.Citruline = subset_samples(ps.Syncom3, Drugs ==  "Citruline")
ps.perc.Citruline2 <- prune_samples(sample_sums(ps.perc.Citruline)>1, ps.perc.Citruline) 

ps.perc.Adrenaline = subset_samples(ps.Syncom3, Drugs ==  "Adrenaline")
ps.perc.Adrenaline2 <- prune_samples(sample_sums(ps.perc.Adrenaline)>1, ps.perc.Adrenaline) 

ps.perc.Pyrone = subset_samples(ps.Syncom3, Drugs ==  "Pyrone")
ps.perc.Pyrone2 <- prune_samples(sample_sums(ps.perc.Pyrone)>1, ps.perc.Pyrone) 

ps.perc.Glutaric = subset_samples(ps.Syncom3, Drugs ==  "Glutaric")
ps.perc.Glutaric2 <- prune_samples(sample_sums(ps.perc.Glutaric)>1, ps.perc.Glutaric) 

ps.perc.Lysine = subset_samples(ps.Syncom3, Drugs ==  "Lysine")
ps.perc.Lysine2 <- prune_samples(sample_sums(ps.perc.Lysine)>1, ps.perc.Lysine) 
#LY9250 LP8539 

#####PCoA for ASV-level data with Bray-Curtis
# Calculates Bray-Curtis distances between samples. Because taxa is in
# columns, it is used to compare different samples. We transpose the
# assay to get taxa to columns
PCoA <- ordinate(ps.perc.Methycitrate2, "PCoA", "bray")
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

### 4 Speceis
#Pyrone no differences
#Methilcitrate Species: 0.002 Complexity 0.001  
#Adrenaline Species: 0.005 Complexity 0.004
#Control Species: 0.001 Complexity 0.001
#Citruline Species: 0.002  Complexity 0.001


#Control: Species: R2=0.34524 F=2.531  P=0.002 
#Complexity  1   0.8729 0.11855 4.3453  0.002 **
#Species     4   1.6693 0.22670 2.0774  0.008 **
#Residual   24   4.8213 0.65476    

#Methy-Citrate Species:   R2 = 0.36381 F=2.6306  P=0.001 
#Complexity  1   0.3934 0.04038 1.4597  0.136    
#Species     4   3.1516 0.32344 2.9233  0.001 ***
#Residual   23   6.1991 0.63619

#Citruline: Species   R2 = 0.37595 F=2.5302  P=0.001
#Complexity  1   0.3980 0.04781 1.6090  0.124    
#Species     4   2.7313 0.32813 2.7605  0.001 ***
#Residual   21   5.1944 0.62405 

#Adrenaline: Species   R2 = 0.32682 F= 2.3303  P= 0.001
#Complexity  1   0.1993 0.02817 1.0043  0.390    
#Species     4   2.1133 0.29865 2.6618  0.001 ***
#Residual   24   4.7636 0.67318 

#Pyrone: Species R2 = 0.34545 F= 1.7944  P=0.009
#Complexity  1   0.3215 0.03943 1.0241  0.006 **   
#Species     4   2.4952 0.30601 1.9869  0.003 **
#Residual   17   5.3372 0.65455 

#Gluraic: Species R2 =  0.2893 1.6282  0.012
#Complexity  1   0.5616 0.14202 2.2936  0.052 .
#Species     2   0.6995 0.17688 1.4283  0.168  
#Residual   11   2.6936 0.68110 

#Lysine: 
#adonis2(formula = ps.3_bray ~ Complexity * Species, data = sample_ps.3c)
#Complexity  1   0.5616 0.14202 2.2936  0.060 .
#Species     2   0.6995 0.17688 1.4283  0.144  
#Residual   11   2.6936 0.68110                

