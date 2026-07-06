packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')
library(ade4)

sapply(packages, require, character.only = TRUE)              

seqtab_Root <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig S16/Fig S16a/seqtab_Root.rds")
seqtab_Frond <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig S16/Fig S16a/seqtab_Frond.rds")

#####PCoA for ASV-level data with Bray-Curtis
PCoA <- ordinate(seqtab_Frond, "PCoA", "bray")
f<-data.frame(PCoA$vectors[, 1:2])
f <- merge(f,seqtab_Frond@sam_data, by=0)
f$cell_layer_n = as.numeric(f$cell_layer_n)

PCoAFrondMicro = aggregate(cbind(Frond.x=Axis.1,Frond.y=Axis.2)~Species*Nutrient,f,mean)
PCoAFrondMicrosd = aggregate(cbind(Frond_ds.x=Axis.1,Frond_ds.y=Axis.2)~Species*Nutrient,f,sd)
Frond <- cbind(PCoAFrondMicro,PCoAFrondMicrosd[,3:4])

PCoA <- ordinate(seqtab_Root, "PCoA", "bray")
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

#Root traits only syncom
x2 = x[x$Microbiome == 'B+',]

#Dimesnsions
nrow(r)
nrow(f)
nrow(x2)

#Conbination
ID = c(166,215,46,234,167,163,37,143,57,149,151,219,165,162,47,153,150,42,225,161,154,
       146,51,218,214,222,236,156,235,58,63,61,142,217,144,221,155,231,141,213,39,152,44,62,
       164,147,168,45,49,59,43,216,36,41)

r2 = cbind(r,ID)
colnames(r2) = c("Row.names","Axis.1","Axis.2","Plate","Column","Row",         
                 "Compartment","Treatment","Nutrient","Species","ID2","cell_layer_n","ID" )
r2x2 = merge(r2, x2[,c(3,5)], by='ID' )
nrow(r2x2)

r2_DIST <- dist(r2x2[31:46,3]) 
r2_DIST_t = as_tibble(r2_DIST, rownames='A')
x2_DIST <- dist(r2x2[31:46,14]) 
x2_DIST_t = as_tibble(x2_DIST, rownames='A')

r2x2_DIST_t = cbind(r2_DIST_t,x2_DIST_t)
colnames(r2x2_DIST_t) = c( 'A','r', 'B','x')

ggplot(data = r2x2_DIST_t, aes(r, x )) +
  geom_point( size = 2,stroke = 1 ) + #, aes(fill = RootComp aes(colour =CCL[,3])
  geom_smooth(method="lm", se=T)+  
  labs(y='Root Microbiome (PCoA1)', x='Root traits (PCoA1)')+
  theme_few()


