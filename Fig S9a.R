packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest",'dada2', "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", "reshape2")
sapply(packages, require, character.only = TRUE)              
setwd("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/")

###Track sequencing abundance
#track <- readRDS("C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/track240-200.rds")
#write.table(track, "track 16S 240-240 Exp1.txt")

# 1. Taxonomy Table
taxa <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/tax_final250_220.rds")

# 2. Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)
#write.table(ps2@sam_data, file = "C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/metadataP.tsv", row.names=FALSE, sep="\t")

# 3. ASV table
asv.table<- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/seqtab_final250_220.rds")
#write.table(asv.table, file = "C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/asv_table.txt", row.names=FALSE, sep="\t")
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)


# 4. Now we can make the phyloseq object
ps2 <- phyloseq(asv.table2, tax_table(taxa), sample_data(meta))

ps.2 = subset_taxa(ps2, Kingdom == "Bacteria" | Kingdom == "Archaea")
ps.pruned <- prune_samples(sample_sums(ps.2)>=2, ps.2)
ps.pruned_taxa <- filter_taxa(ps.2, function(x) sum(x) > .005, TRUE)

ps.3 = subset_samples(ps.pruned, Compartment !=  "BLANK")
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))

#  Plotting Relative Abundance Bar Charts####
# phylum-level
ps.perc <- transform_sample_counts(ps.3, function(x) x / sum(x)) 
ps.perc.Front <- transform_sample_counts(ps.3.Front, function(x) x / sum(x)) 
ps.perc.Root <- transform_sample_counts(ps.3.Root.H4, function(x) x / sum(x)) 
ps.perc.Water <- transform_sample_counts(ps.3.Water, function(x) x / sum(x)) 

taxa_level <- function(physeq,which_level){
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  OTU <- otu_table(physeq)
  SAM <- sample_data(physeq)
  OTU_taxonomy <- tax_table(physeq)
  new_abund_table<-NULL
  if(which_level=="Otus"){
    OTU_tree <- phy_tree(physeq)
    new_abund_table<-OTU
  } else {
    list<-na.omit(unique(OTU_taxonomy[,which_level]))
    new_abund_table<-NULL
    for(i in list){
      rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
      tmp<-data.frame(rowSums(OTU[,rt]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  OTU<-as.data.frame(as(new_abund_table,"matrix"))
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(SAM)
  #reconstruct the phyloseq object
  physeq<-NULL
  if(which_level=="Otus"){
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }
  return(physeq)
}

ps.phyla.perc <- taxa_level(ps.perc, "Phylum")
ps.phyla.perc.Front <- taxa_level(ps.perc.Front, "Phylum")
ps.phyla.perc.Root <- taxa_level(ps.perc.Root, "Phylum")
ps.phyla.perc.Water <- taxa_level(ps.perc.Water, "Phylum")

# identify the 10 most abundant phylum
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:12])
phylum.10.Front <- names(sort(taxa_sums(ps.phyla.perc.Front), TRUE)[1:12])
phylum.10.Root <- names(sort(taxa_sums(ps.phyla.perc.Root), TRUE)[1:12])
phylum.10.Water <- names(sort(taxa_sums(ps.phyla.perc.Water), TRUE)[1:12])

# now we can use the list of the top phylum and subset those from the phyloseq object
ps.phylum.10 <- prune_taxa(phylum.10, ps.phyla.perc)
ps.phylum.10.Front <- prune_taxa(phylum.10.Front, ps.phyla.perc.Front)
ps.phylum.10.Root <- prune_taxa(phylum.10.Root, ps.phyla.perc.Root)
ps.phylum.10.Water <- prune_taxa(phylum.10.Water, ps.phyla.perc.Water)

melt.phylum <- psmelt(ps.phylum.10)
melt.phylum.Front <- psmelt(ps.phylum.10.Front)
melt.phylum.Root <- psmelt(ps.phylum.10.Root)
melt.phylum.Water <- psmelt(ps.phylum.10.Water)

safe_colorblind_palette <- rev(c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                 "#44AA99", "#999933", "#882255", "#661100", "#CC6677", "#DDCC77"))

#cargar dplyr
library(dplyr)

#Calcular medias por grupos de todas las columnas
melt.phylum.mean <- melt.phylum%>%group_by(Compartment, OTU)%>%
  summarise_all(mean)

melt.phylum.mean.Root <- melt.phylum.Root%>%group_by(Specie.1, OTU)%>%
  summarise_all(mean)

melt.phylum.mean.Front <- melt.phylum.Front%>%group_by(Specie.1, OTU)%>%
  summarise_all(mean)

melt.phylum.mean.Water <- melt.phylum.Water%>%group_by(Specie.1, OTU)%>%
  summarise_all(mean)


###Graficos
melt.phylum.mean.Root$Specie.1 <- factor(melt.phylum.mean.Root$Specie.1,rev(c('PS','LY5280','SI9227','SP7498','SP9509',
                                                                              'SI7820','SP5543','SP9192','LM7200','LT9243',
                                                                              'LP7760','LV7650','LJ9250','LP0049','LM8389','LT9109',
                                                                              'LA7339','LM7123','LP8539', 'Control', 'Inoculum')))
melt.phylum.mean.Front$Specie.1 <- factor(melt.phylum.mean.Front$Specie.1,rev(c('PS','LY5280','SI9227','SP7498','SP9509',
                                                                                'SI7820','SP5543','SP9192','LM7200','LT9243',
                                                                                'LP7760','LV7650','LJ9250','LP0049','LM8389','LT9109',
                                                                                'LA7339','LM7123','LP8539')))
melt.phylum.mean.Water$Specie.1 <- factor(melt.phylum.mean.Water$Specie.1,rev(c('PS','LY5280','SI9227','SP7498','SP9509',
                                                                                'SI7820','SP5543','SP9192','LM7200','LT9243',
                                                                                'LP7760','LV7650','LJ9250','LP0049','LM8389','LT9109',
                                                                                'LA7339','LM7123','LP8539', 
                                                                                'Control', 'Inoculum')))

melt.phylum.mean$Compartment   <- factor(melt.phylum.mean$Compartment, c('Inoculum','Front', 'Root','Water'))
melt.phylum.mean.Water$Compartment   <- factor(melt.phylum.mean.Water$Compartment, c('Inoculum','Front', 'Root','Water'))
melt.phylum.mean.Front$Compartment   <- factor(melt.phylum.mean.Front$Compartment, c('Inoculum','Front', 'Root','Water'))
melt.phylum.mean.Root$Compartment   <- factor(melt.phylum.mean.Root$Compartment, c('Inoculum','Front', 'Root','Water'))

melt.phylum.mean$OTU <- factor(melt.phylum.mean$OTU,c("Proteobacteria","Cyanobacteria","Bacteroidetes", "Firmicutes","Actinobacteria","Gemmatimonadetes","Acidobacteria","TM6_(Dependentiae)","Chlamydiae","Verrucomicrobia","Spirochaetae","Planctomycetes"))
melt.phylum.mean.Front$OTU <- factor(melt.phylum.mean.Front$OTU,rev(c("Proteobacteria","Cyanobacteria","Bacteroidetes","Actinobacteria","Acidobacteria","Firmicutes","Gemmatimonadetes","TM6_(Dependentiae)","Verrucomicrobia","Planctomycetes","Parcubacteria","Chloroflexi" )))
melt.phylum.mean.Root$OTU <- factor(melt.phylum.mean.Root$OTU,rev(c("Proteobacteria","Cyanobacteria", "Bacteroidetes",  "Actinobacteria","Gemmatimonadetes", "Acidobacteria", "Planctomycetes", "Spirochaetae", "TM6_(Dependentiae)", "Firmicutes","Chlamydiae","Verrucomicrobia")))
melt.phylum.mean.Water$OTU <- factor(melt.phylum.mean.Water$OTU,rev(c("Proteobacteria","Cyanobacteria","Bacteroidetes","Actinobacteria","Gemmatimonadetes","Verrucomicrobia","Firmicutes","Acidobacteria","Chlamydiae","Planctomycetes","TM6_(Dependentiae)","Spirochaetae")))


ggplot(melt.phylum.mean, aes(y = Compartment , x = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(melt.phylum.mean.Root, aes(y = Specie.1 , x = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + labs(title ="Root")+theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(melt.phylum.mean.Front, aes(y = Specie.1 , x = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + labs(title ="Front")+theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(melt.phylum.mean.Water[melt.phylum.mean.Water$Specie.1!='Control',], aes(y = Specie.1 , x = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + labs(title ="Water", legend ="Phylum")+theme(axis.text.x = element_text(angle = 45, hjust=1))

