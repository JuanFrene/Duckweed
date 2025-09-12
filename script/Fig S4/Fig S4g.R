packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest",'dada2', "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", "reshape2", 'heatmaply','ggcorrplot',
              "viridis",'dendextend','dplyr','cowplot','ggdendro','pals','paletteer')
sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/Exp1/")

TablaCompleta = read.table("Root traits.txt", header = TRUE, row.names = 1)

TablaCompleta$Subsample = as.factor(TablaCompleta$Subsample)

Root.pca <- prcomp(TablaCompleta[,-c(1:4)], scale = TRUE)
x<-data.frame(Root.pca$x[,1:2])
x <- cbind(TablaCompleta[,1:4],x )
rownames(x)<-x$Row.names
x$Row.names <- NULL

PCA_Root_mean = aggregate(cbind(mean.x=PC1,mean.y=PC2)~Species*Microbiome ,x,mean)
PCA_Root_sd = aggregate(cbind(ds.x=PC1,ds.y=PC2)~Species*Microbiome ,x,sd)
Root_trait2 <- cbind(PCA_Root_mean,PCA_Root_sd[,-(1:2)])

setwd("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/")

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

ps.3.Frond = subset_samples(ps.3, Compartment ==  "Frond")



#####PCoA for ASV-level data with Bray-Curtis
PCoA.F <- ordinate(ps.3.Frond, "PCoA", "bray")
f<-data.frame(PCoA.F$vectors[, 1:2])
f <- merge(f,ps.3.Frond@sam_data, by=0)
rownames(f)<-f$Row.names
f$Row.names <- NULL
nrow(f)

PCoA = aggregate(cbind(Frond.x=Axis.1,Frond.y=Axis.2)~Specie.1,f,mean)
PCoAsd = aggregate(cbind(Frond.x=Axis.1,Frond.y=Axis.2)~Specie.1,f,sd)
Frond <- merge(PCoA,PCoAsd, by="Specie.1")

all = merge(Root_trait2[Root_trait2$Microbiome=='B+',],Frond, by='Species' )


ggplot(data = all, aes(mean.x, Frond )) +
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax =mean.x + ds.x ),size = 0.2, alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =Root.1 - Root.sd.1,ymax = Root.1 + Root.sd.1),size = 0.2, alpha = 0.4)  + 
  geom_point( size = 3,stroke = 1) + #, aes(fill = RootComp
  geom_smooth(method="lm", se=T)+  
  labs(y='Root Microbiome (PCoA1)', x='Root traits (PCoA1)')+
  scale_fill_gradientn(colours=c('blue','red'))+
  theme_few()+
  xlim(NA,3)

xB = x[x$Microbiome=='B+',]
Plant = c(106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200)

charla.F = merge(xB2,f, by ='Plant')
root.abund = vegdist(charla.F$PC1, method = "euclidean")
frondmicro.abund2 = vegdist(charla.F$Axis.2, method = "bray")
Root.frond2 = mantel(root.abund, frondmicro.abund2, method = "spearman", permutations = 9999, na.rm = TRUE)
#Mantel statistic r: -0.05472  Significance: 0.8791
