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
Alfadiv2$Shannon = as.numeric(Alfadiv2$Shannon)


Alfadiv2$Root = factor(Alfadiv2$Root, c('High', 'Medium', 'Low'))
Alfadiv2$Species = factor(Alfadiv2$Species, c('Control','LP8539','LA7339','LP0049' ,'LJ9250' ,'LT9243','SP9509','SI7820','SP7498','SI9227'))

ggplot(Alfadiv2, aes(Root, Shannon), ) + 
  geom_boxplot(aes(fill=Complexity)) +
  facet_grid(.~Species)+
  theme_few()+
  scale_fill_gradientn(colours=c('blue','red'))+ 
  labs(x='Synthetic Root', y = 'Shannon alphadiversity') + 
  theme(legend.position="top",axis.text.x = element_text(angle=65, vjust=0.6))



summary(aov(Shannon ~ Root  + (1|Sample), data = Alfadiv2[Alfadiv2$Specie=='LA7339',]))


shannon.model.Root<-lmer(Shannon ~ Root*Specie  + (1|Sample), data = Alfadiv2)
anova(shannon.model.Root)
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Root        0.86897 0.43449     2 84.246 28.4446 3.638e-10 ***
#Specie      1.30756 0.14528     9 84.197  9.5114 7.368e-10 ***
#Root:Specie 0.59946 0.03330    18 84.160  2.1803  0.009126 **

#Control 
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Root 0.46576 0.23288     2 5.8815  28.288 0.0009605

#LA7339 
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.1111 0.055552     2 5.2318  3.6313 0.1026

#LJ9250
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.01205 0.0060248     2 2.6004  0.5969 0.6119

#LP0049 
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.13909 0.069544     2 6.3235  8.6076 0.01567 *

#LP8539
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.16122 0.080612     2 6.7364  11.117 0.007348 **

#LT9243
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.058018 0.029009     2 5.8109  1.2865 0.3447

#SI7820 
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.11087 0.055434     2 3.3858  4.8604 0.1011


#SI9227
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.089519 0.044759     2 4.942  7.1779 0.03453 *

#SP7498
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.0079697 0.0039848     2 6.577  0.2209 0.8075


#SP9509
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Root 0.020626 0.010313     2 6.5431  0.7248 0.5196

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

#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF DenDF  F value    Pr(>F)    
#Root     0.99539 0.49769     2 106.7  19.313    6.941e-08#
  
  
  