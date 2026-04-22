###Heatmap
library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggplot2)
library(ggthemes)
library(ggdendro)
library(pals)
library(paletteer)
require(GGally)
require(CCA)
require(CCP)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 1/Fig 1e/")

Table = read.table("Root traits.txt", header = TRUE, row.names = 1)

Table$Subsample = as.factor(Table$Subsample)

varespec.bray <- vegdist(Table[,-c(1:4)], method = "euclidean") # dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan
pcoaVS <- pco(varespec.bray, negvals = "zero", dround = 0)
pcoaVS$vectors[,1:2]

f<-data.frame(pcoaVS$vectors[,1:2])
f <- cbind(Table[,c(1:3,23)],f )
rownames(f)<-x$Row.names
f$Row.names <- NULL

PCA_Comp_mean = aggregate(cbind(mean.x=X1,mean.y=X2,mean.z=CCL)~Species*Microbiome ,f[-c(190),],mean)
PCA_Comp_meands = aggregate(cbind(ds.x=X1,ds.y=X2)~Species*Microbiome,f[-(190),],sd)
PCoA_Comp_mean2 =cbind(PCA_Comp_mean,PCA_Comp_meands[,c(3:4)])


paleta_alive <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                  "#14AA99","#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                  "#49AA91", "#999983", "#882249", "#661830", "#6699DC", "#828389")

ggplot(data = PCoA_Comp_mean2,aes(mean.x,mean.y)) + #PCoA_Comp.F_mean2 mean.x mean.y
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),linewidth = 0.1,alpha = 0.3) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),linewidth = 0.01,alpha = 0.1)  + 
  geom_point(size=2.5, aes(color = mean.z,  shape=Microbiome),stroke = 1) + #size = 4,
  labs(x='PCoA1', y='PCoA2') +
  scale_color_gradientn(colours=c('blue','red'))+
  theme_few()



