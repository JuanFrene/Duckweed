packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'reshape2')
sapply(packages, require, character.only = TRUE)              
library(dada2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Microtubulos/16S/")

###ASV table
asv.table2 <- readRDS('G:/My Drive/labs/Nottingham/Duckweed/Exp Microtubulos/16S/seqtab_final.rds')
asv.table<- otu_table(asv.table2, taxa_are_rows=FALSE)

##ASing to the syncom
taxa <- assignTaxonomy(asv.table, "G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/Syncom_Sequence.pcr.unique.fasta", multithread=TRUE)#, minBoot =97

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

#Now we can make the phyloseq object
ps <- phyloseq(asv.table, tax_table(taxa), sample_data(meta))
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps.1 = subset_taxa(ps, Genus  != "NA")

ps.2 = subset_samples(ps.1, ID !=  "Blanck")
ps.2.pruned <- prune_taxa(taxa_sums(ps.2)>=1, ps.2)
ps.2.pruned2 <- prune_samples(sample_sums(ps.2.pruned)>0, ps.2.pruned) 

ps.3 = subset_samples(ps.2, Root !=  "OR")
ps.3.pruned <- prune_taxa(taxa_sums(ps.3)>=1, ps.3)
ps.3.pruned2 <- prune_samples(sample_sums(ps.3.pruned)>0, ps.3.pruned) 

ps.4 = subset_samples(ps.3, ID !=  "Syncom")
ps.4.pruned <- prune_taxa(taxa_sums(ps.4)>=1, ps.4)
ps.4.pruned2 <- prune_samples(sample_sums(ps.4.pruned)>0, ps.4.pruned) 

ps.SR2 = subset_samples(ps.4.pruned2, Compartment ==  "SR")
ps.SR1 <- prune_taxa(taxa_sums(ps.SR2)>=1, ps.SR2)
ps.SR <- prune_samples(sample_sums(ps.SR1)>0, ps.SR1) 
ps.SR.H <- transform(ps.SR, "hellinger")
ps.SR.perc <- transform_sample_counts(ps.SR, function(x) x / sum(x)) 

####Alfadiversity
ordered(sample_sums(ps.SR))

# let's calcuolate Shannon diversity and Richness (Observed)
alpha.div <- estimate_richness(ps.SR, measures=c("Shannon", "Observed",'Chao', 'Simpson'))
even <- evenness(ps.SR, 'pielou')

Alfadiv <- ps.SR@sam_data

# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses
Alfadiv$Shannon <- paste(alpha.div$Shannon)
Alfadiv$Specie <- as.factor(Alfadiv$Species)
Alfadiv$Root <- as.factor(Alfadiv$Root)
Alfadiv$Species <- as.factor(Alfadiv$Species)
Alfadiv$Sample <- as.numeric(Alfadiv$Sample)

Alfadiv = data.frame(Alfadiv)

Alfadiv2 = Alfadiv[Alfadiv$Shannon > 0.01,]

shannon.model.Root<-lmer(Shannon ~ Specie  + (1|Sample), data = Alfadiv2[Alfadiv2$Root=='Low',])
aov = anova(shannon.model.Root)

#High
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Specie 0.35904 0.039893     9 23.555  2.5697 0.03196 *

#Medium
#Type III Analysis of Variance Table with Satterthwaite's method
#       Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)   
#Specie 0.57581 0.063979     9 24.568  4.5023 0.00142

#Low
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Specie 0.91438  0.1016     9 29.131  5.3306 0.0002534 ***
---
  
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF DenDF  F value    Pr(>F)    
#Root     0.99539 0.49769     2 106.7  19.313    6.941e-08
  
  