packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')
sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/")

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)
meta3[1:6,1:5]

###ASV table
asv.table <- readRDS('G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/seqtab_final240_220.rds')
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

##ASing to the syncom
taxaSyncom <- assignTaxonomy(asv.table2, "G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/Syncom_Sequence.pcr.unique.fasta", multithread=TRUE)#, minBoot = 98

#Now we can make the phyloseq object
ps.Syncom <- phyloseq(asv.table2, tax_table(taxaSyncom), sample_data(meta))
taxa_names(ps.Syncom) <- paste0("ASV", seq(ntaxa(ps.Syncom)))
ps.Syncom = subset_taxa(ps.Syncom, Genus  != "NA")

ps.Syncom.3 = subset_samples(ps.Syncom, Treatment !=  "ConNeg")
ps.pruned <- prune_taxa(taxa_sums(ps.Syncom.3)>=1, ps.Syncom.3)

ps.Syncom.3.H <- transform(ps.pruned, "hellinger")
ps.Syncom.3.perc <- transform_sample_counts(ps.pruned, function(x) x / sum(x)) 

ps.Syncom.SC = subset_samples(ps.Syncom.3.H, Species !=  "Control")
ps.Syncom.SC2 = subset_samples(ps.Syncom.SC, Species !=  "Inoculum")

ps.Syncom.NOSyncom = subset_samples(ps.Syncom.SC2, Treatment !=  "SynCom")
ps.Syncom.NOSyncom2 <- prune_taxa(taxa_sums(ps.Syncom.NOSyncom)>0, ps.Syncom.NOSyncom)

ps.Syncom.Syncom = subset_samples(ps.Syncom.SC2, Treatment ==  "SynCom")
ps.Syncom.Syncom2 <- prune_taxa(taxa_sums(ps.Syncom.Syncom)>1, ps.Syncom.Syncom)

ps.Syncom.Syncom2 <- prune_taxa(taxa_sums(ps.Syncom.Syncom)>1, ps.Syncom.Syncom)
ps.Syncom.Syncom3 <- prune_samples(sample_sums(ps.Syncom.Syncom2)>1, ps.Syncom.Syncom2) 

ps.Syncom.Front = subset_samples(ps.Syncom.Syncom, Compartment ==  "Front")
ps.Syncom.Root = subset_samples(ps.Syncom.Syncom, Compartment ==  "Root")

ps.Syncom.Front2 <- prune_taxa(taxa_sums(ps.Syncom.Front)>0, ps.Syncom.Front)
ps.Syncom.Front3 <- prune_samples(sample_sums(ps.Syncom.Front)>0, ps.Syncom.Front) 
ps.Syncom.Front4 <- subset_samples(ps.Syncom.Front3, Species !=  "PS") 

ps.Syncom.Root2 <- prune_taxa(taxa_sums(ps.Syncom.Root)>0, ps.Syncom.Root)
ps.Syncom.Root3 <- prune_samples(sample_sums(ps.Syncom.Root2)>0, ps.Syncom.Root2) 
ps.Syncom.Root4 <- subset_samples(ps.Syncom.Root3, Species !=  "PS") 

ps.Syncom.High = subset_samples(ps.Syncom.Root4, Nutrient ==  "High")
ps.Syncom.Medium = subset_samples(ps.Syncom.Root4, Nutrient ==  "Medium")
ps.Syncom.Low = subset_samples(ps.Syncom.Root4, Nutrient ==  "Low")

#####PCoA for ASV-level data with Bray-Curtis
PCoA <- ordinate(ps.Syncom.Front4, "PCoA", "bray")
f<-data.frame(PCoA$vectors[, 1:2])
f <- merge(f,ps.Syncom.Front3@sam_data, by=0)
f$cell_layer_n = as.numeric(f$cell_layer_n)

PCoAFrondMicro = aggregate(cbind(Frond.x=Axis.1,Frond.y=Axis.2)~Species*Nutrient,f,mean)
PCoAFrondMicrosd = aggregate(cbind(Frond_ds.x=Axis.1,Frond_ds.y=Axis.2)~Species*Nutrient,f,sd)
Frond <- cbind(PCoAFrondMicro,PCoAFrondMicrosd[,3:4])

PCoA <- ordinate(ps.Syncom.Root4, "PCoA", "bray")
r<-data.frame(PCoA$vectors[, 1:2])
r <- merge(r,ps.Syncom.Root3@sam_data, by=0)
r$cell_layer_n = as.numeric(r$cell_layer_n)

PCoARootMicro = aggregate(cbind(Root.x=Axis.1,Root.y=Axis.2)~Species*Nutrient,r,mean)
PCoARootMicrosd = aggregate(cbind(Root_ds.x=Axis.1,Root_ds.y=Axis.2)~Species*Nutrient,r,sd)
Root <- cbind(PCoARootMicro,PCoARootMicrosd[,3:4])

Micro <- cbind(Root,Frond[,-c(1:2)])

setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/Exp2/")

TablaCompleta = read.table("Root_trait_exp2.txt", header = TRUE, row.names = 1)
TablaCompleta[1:5,1:5]

Root.pca <- prcomp(TablaCompleta[,6:25], scale = TRUE)
x<-data.frame(Root.pca$x[,1:2])
x <- cbind(TablaCompleta[,1:4],x )
rownames(x)<-x$Row.names
x$Row.names <- NULL

PCA_Root_mean = aggregate(cbind(mean.x=PC1,mean.y=PC2)~Species*Nutrition ,x,mean)
PCA_Root_sd = aggregate(cbind(ds.x=PC1,ds.y=PC2)~Species*Nutrition ,x,sd)
Root_trait2 <- cbind(PCA_Root_mean,PCA_Root_sd[,-(1:2)])

All <- cbind(Micro,Root_trait2[,-c(1:2)])

CCL = aggregate(cell_layer_n~Species*Nutrition ,TablaCompleta,mean)
All2 <- cbind(All,CCL[,3])

All$Nutrient = factor(All$Nutrient, c('High', 'Medium', 'Low'))
ggplot(data = All2, aes(mean.x, Root.x )) +
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax =mean.x + ds.x ),size = 0.2, alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =Root.x - Root_ds.x,ymax = Root.x + Root_ds.x),size = 0.2, alpha = 0.4)  + 
  geom_point( size = 3,stroke = 1, aes(colour =CCL[,3])) + #, aes(fill = RootComp
  facet_grid(.~Nutrient, scales = "fixed")+
  geom_smooth(method="lm", se=T)+  
  labs(y='Root Microbiome (PCoA1)', x='Root traits (PCoA1)')+
  scale_color_gradientn(colours=c('blue','red'))+
  theme_few()
  
ggplot(data = All2, aes(mean.x, Frond.x )) +
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax =mean.x + ds.x ),size = 0.2, alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =Frond.x - Frond_ds.x,ymax = Frond.x + Frond_ds.x),size = 0.2, alpha = 0.4)  + 
  geom_point( size = 3,stroke = 1, aes(colour =CCL[,3])) + #, aes(fill = RootComp
  facet_grid(.~Nutrient,scales = "free")+
  geom_smooth(method="lm", se=T)+  
  labs(y='Frond Microbiome (PCoA1)', x='Root traits (PCoA1)')+
  scale_color_gradientn(colours=c('blue','red'))+
  theme_few()

Front.micro4 = vegdist(log(All[All$Nutrient=='Low',]$Frond.x+2), method = "bray")
Root.micro4 = vegdist(log(All[All$Nutrient=='Low',]$Root.x+2), method = "bray")
root.trait4 = vegdist(All[All$Nutrient=='Low',]$mean.x, method = "euclidea")
  
Root.root6 = mantel(Root.micro6, root.trait6, method = "spearman", permutations = 9999, na.rm = TRUE)
Frond.root6 = mantel(Front.micro6, root.trait6, method = "spearman", permutations = 9999, na.rm = TRUE)

Root.root5 = mantel(Root.micro5, root.trait5, method = "spearman", permutations = 9999, na.rm = TRUE)
Frond.root5 = mantel(Front.micro5, root.trait5, method = "spearman", permutations = 9999, na.rm = TRUE)

Root.root4 = mantel(Root.micro4, root.trait4, method = "spearman", permutations = 9999, na.rm = TRUE)
Frond.root4 = mantel(Front.micro4, root.trait4, method = "spearman", permutations = 9999, na.rm = TRUE)


#Low Root:r: 0.3643 P: 0.091667 , frond:r:-0.1857, P = 0.70417 
#Medium Root:r: -0.1321  P:0.65694  , frond:r:0.1464, P = 0.21944  
#High Root:r: 0.025 P: 0.44306  , frond:r:0.01786 , P = 0.40417  
