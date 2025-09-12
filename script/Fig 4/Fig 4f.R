packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')

sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

# Control sequence ####
#track <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/track.rds")
#write.table(track, "track 16S.txt")


###ASV table
asv.table <- readRDS('G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/seqtab_final.rds')
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

##ASing to the syncom
taxaSyncom <- assignTaxonomy(asv.table2, "G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/Syncom_Sequence.pcr.unique.fasta", multithread=TRUE)#, minBoot = 98

#Metadata
meta <- read.table("Metadata SR.txt", header = TRUE, row.names = 1)

#Now we can make the phyloseq object
ps.Syncom <- phyloseq(asv.table2, tax_table(taxaSyncom), sample_data(meta))

taxa_names(ps.Syncom) <- paste0("ASV", seq(ntaxa(ps.Syncom)))
ps.Syncom = subset_taxa(ps.Syncom, Genus  != "NA")

ps.pruned <- prune_taxa(taxa_sums(ps.Syncom)>=1, ps.Syncom)
ps.pruned3 <- prune_samples(sample_sums(ps.pruned)>1, ps.pruned) 
ps.perc <- transform_sample_counts(ps.pruned3, function(x) x / sum(x)) 

ps.perc.Control = subset_samples(ps.perc, Drugs ==  "No")
ps.perc.Control2 <- prune_samples(sample_sums(ps.perc.Control)>0, ps.perc.Control) 


ps.perc.Lysine = subset_samples(ps.perc, Drugs ==  "Lysine")
ps.perc.Lysine2 <- prune_samples(sample_sums(ps.perc.Lysine)>0, ps.perc.Lysine) 

#####PCoA for ASV-level data with Bray-Curtis
PCoA <- ordinate(ps.perc.Control2, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x, ps.perc@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Species*Root,x,mean)
PCoAmeands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Species*Root,x,se)
PCoAmeanLJ =merge(PCoAmean[c(1,4),],PCoAmeands[c(1,4),-(1)],by="Root")
PCoAmeanSP =merge(PCoAmean[-c(1,4),],PCoAmeands[-c(1,4),-(1)],by="Root")
PCoAmean2 = rbind(PCoAmeanLJ,PCoAmeanSP)


paleta_alive <- c("#C200FF",'#8B7500','#FF0000','#00008B',"#FFB919",'#FF7F50',"#00CC7A")

PCoAmean2$Root = factor(PCoAmean2$Root, c('High', 'Medium','Low'))
ggplot(data = PCoAmean2, aes(mean.x, mean.y, label=Species)) +
    geom_vline(xintercept = 0, size = 2, color = "#D9D9D9",linetype = "longdash")+
    geom_hline(yintercept = 0, size = 2, color = "#D9D9D9",linetype = "longdash")+
    geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),linewidth = 0.1, alpha = 0.15) +
    geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),linewidth = 0.1, alpha = 0.15)  + 
    geom_point( size = 6, aes(colour = Root, shape = Species),stroke = 1) + 
    #scale_fill_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
    labs(y=ylab_text, x=xlab_text)+ theme_few()+
    scale_colour_manual(values = paleta_alive)

ggplot(x, aes(Axis.1, Axis.2)) + #[PCoAmean$Specie!='PS',]
  geom_point(aes(color=Species,shape=Root, size=6)) +
  #geom_text(hjust=0.5, vjust=-0.5) + 
  theme_few() +
  theme(legend.position="top")+ 
  scale_fill_paletteer_d("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  #scale_colour_gradientn(colours=c("Blue","Red"))+
  labs(title = "Bacterial community associated with Duckweed Root",
       y=ylab_text, 
       x=xlab_text)    

library(paletteer)
library(pals)


  ###PERMANOVA
# Calculate bray curtis distance matrix
ps.3_bray <- phyloseq::distance(ps.perc.Control2, method = "bray")
# make a data frame from the sample_data
sample_ps.3c <- data.frame(sample_data(ps.perc.Control))
# Adonis test
adonis2(ps.3_bray ~ Root*Species , data = sample_ps.3c)

#Control            Df SumOfSqs      R2      F Pr(>F)  
#Root          2   1.1450 0.25706 2.6632  0.020 *
#Species       1   0.1273 0.02858 0.5923  0.677  
#Root:Species  1   0.1724 0.03870 0.8019  0.519  
#Residual     14   3.0095 0.67565                
#Total        18   4.4541 1.00000 


ps.Lysine_bray <- phyloseq::distance(ps.perc.Lysine2, method = "bray")
# make a data frame from the sample_data
sample_ps.Lysine <- data.frame(sample_data(ps.perc.Lysine2))
# Adonis test
adonis2(ps.Lysine_bray ~ Root*Species , data = sample_ps.Lysine)

#Lysine: 
#            Df SumOfSqs      R2      F Pr(>F)   
#Root          2   0.4731 0.08792 1.0734  0.359   
#Species       1   0.8163 0.15170 3.7041  0.003 **
#Root:Species  2   0.3453 0.06416 0.7834  0.662   
#Residual     17   3.7465 0.69622                 
#Total        22   5.3813 1.00000   

