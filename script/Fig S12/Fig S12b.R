packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2','dplyr')


setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/")

###Metadata
seqtabExp2 <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig S12/Fig S12b/seqtab Exp2.rds")

ps.Syncom.SC = subset_samples(seqtabExp2, Species !=  "Control")
ps.Syncom.Syncom = subset_samples(ps.Syncom.SC, Treatment !=  "NB")
ps.Syncom.Syncom2 = subset_samples(ps.Syncom.Syncom, Treatment !=  "NatCom")

#Unique taxa
Syncom_taxa = data.frame(ps.Syncom.Syncom2@tax_table)
unique(Syncom_taxa$Order)
##Genus 25
##family 19 
#Order 14
##Class 6 
##Phylum 4

#  Plotting Relative Abundance Bar Charts####
# phylum-level
ps.perc <- transform_sample_counts(ps.3, function(x) x / sum(x)) 
ps.phyla.perc <- taxa_level(ps.perc, "Phylum")
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:12])
ps.phylum.10 <- prune_taxa(phylum.10, ps.phyla.perc)
melt.phylum <- psmelt(ps.phylum.10)


############Natural community
seqtabExp1 <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig S12/Fig S12b/seqtab Exp1.rds")

ps.3.Root = subset_samples(seqtabExp1, Compartment ==  "Root")
ps.3.Root.pruned <- prune_taxa(taxa_sums(ps.3.Root)>=5, ps.3.Root)


#Unique taxa
NatCom_taxa = data.frame(ps.3.Root.pruned@tax_table)
unique(NatCom_taxa$Family)
##Genus 387/R=303
##family 192/R= 159
#Order 103/R=85
##Class 59 /R =49
##Phylum 26/R=22


#  Plotting Relative Abundance Bar Charts####
# phylum-level
ps.perc.natcom <- transform_sample_counts(ps.3.SCSI.pruned, function(x) x / sum(x)) 
ps.phyla.perc.natcom <- taxa_level(ps.perc.natcom, "Order")
phylum.10.natcom <- names(sort(taxa_sums(ps.phyla.perc.natcom), TRUE))
ps.phylum.10.natcom <- prune_taxa(phylum.10.natcom, ps.phyla.perc.natcom)
melt.phylum.natcom <- psmelt(ps.phylum.10.natcom)
melt.phylum.mean.natcom <- data.frame(melt.phylum.natcom%>%group_by(Compartment, OTU)%>%summarise_all(mean))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#CC6677", "#DDCC77", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#CC6677", "#DDCC77")

safe_colorblind_palette2 <- c("black", "black", "black", "black", "black", "black", 
                             "black", "black", "black", "black", "black", "black")

Phylum_order = c("Proteobacteria", "Bacteroidetes","Actinobacteria","Firmicutes","Gemmatimonadetes","Cyanobacteria", "Acidobacteria", "Planctomycetes", "Chlorobi", "TM6_(Dependentiae)", "Chlamydiae","Verrucomicrobia")
Family_order = c('Alcaligenaceae','Bacillaceae','Burkholderiaceae','Caulobacteraceae','Comamonadaceae',
                'Dermacoccaceae','Micrococcaceae','Nocardiaceae','Paenibacillaceae','Planococcaceae','Pseudomonadaceae',
                'Rhizobiaceae','Sphingomonadaceae','Spirosomaceae','Staphylococcaceae','Streptomycetaceae','Xanthobacteraceae','Yersiniaceae')
Order_order = c("Corynebacteriales", "Sphingomonadales","Hyphomicrobiales",
                "Rhizobiales","Micrococcales","Caulobacterales",  
                "Bacillales","Paenibacillales","Streptomycetales", 
                "Enterobacterales","Staphylococcales","Burkholderiales",  
                "Pseudomonadales","Cytophagales"  )


melt.phylum.mean.natcom$OTU <- factor(melt.phylum.mean.natcom$OTU,rev(Order_order))

melt.phylum.mean233 = na.omit(melt.phylum.mean.natcom[,c(1,2,4)])

melt.phylum.mean233$Compartment = factor(melt.phylum.mean233$Compartment, c('Water','Root', 'Front'))
ggplot(melt.phylum.mean233, aes(x = Compartment , y = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette2) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(0,0.5) + theme(legend.position = "none")
  
Syncom = data.frame(t(c('Syncom', 0.945)))
colnames(Syncom) = c('Trea','Syncom')
Syncom$Syncom = as.numeric(Syncom$Syncom)

ggplot(Syncom, aes(x = Trea, y = Syncom )) + theme_few() +
  geom_bar(stat = "identity")+
  ylim(0,1) + ylab('Abundenace')+ xlab('')+
  theme(legend.position = "none")

