###PCoA Exudates
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

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exudados/")

Exudados <- read.table("Exudados.txt", header = TRUE, row.names = 1)
Exudados[,1:5]

varespec.bray <- vegdist(Exudados[,-c(1:4)], method = "bray") # -c(3,27)dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan
pcoaVS <- pcoa(varespec.bray)
pcoaVS$vectors[,1:2]

f<-data.frame(pcoaVS$vectors[,1:4])
f <- cbind(Exudados[,c(1:4)],f )#-c(3,27)
rownames(f)<-x$Row.names
f$Row.names <- NULL

PCA_Comp_mean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2,CCL=CCL,CCL_ds=CCL_ds)~Species ,f,mean)
PCA_Comp_meands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Species,f,sd)
PCoA_Comp_mean2 =cbind(PCA_Comp_mean,PCA_Comp_meands[,c(2:3)])

ggplot(data = PCoA_Comp_mean2, aes(CCL, mean.x, )) +
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.2, alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =mean.x - ds.x ,ymax = mean.x + ds.x),size = 0.2, alpha = 0.4)  + 
  geom_hline(yintercept = 0,size = 0.2,color = "black",linetype = "longdash")+
  geom_point( size = 2, aes(colour = Species),stroke = 1) + 
  geom_smooth(method="lm",size = 0.5, se=T,,color='black') + 
  labs(y='PCoA 1', x='CCL') +
  theme_few()+
  scale_color_manual(values = paleta_alive)
