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

######Set work directory
setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

TablaCompleta = read.table("Root trait LC50.txt", header = TRUE, row.names = 1)

TablaCompleta2 = TablaCompleta[TablaCompleta$Treatment!='C5',]

varespec.bray <- vegdist(TablaCompleta2[,-c(1:3,5)], method = "bray") # dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan
pcoa <- pcoa(varespec.bray)
pcoa$vectors[,1:2]

f<-data.frame(pcoa$vectors[,1:2])
f <- cbind(TablaCompleta2[,c(1:3)],f )
rownames(f)<-f$Row.names
f$Row.names <- NULL

PCA_Comp_mean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Treatment*Dosis, f[-7,],mean)
PCA_Comp_meands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Treatment*Dosis,f[-7,],se)
PCoA_Comp_mean2 =cbind(PCA_Comp_mean,PCA_Comp_meands[,c(3:4)])


paleta_alive <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                  "#14AA99","#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                  "#49AA91", "#999983", "#882249", "#661830", "#6699DC", "#828389")




ggplot(data = PCoA_Comp_mean2,aes(mean.x,mean.y)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),linewidth = 0.1,alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),linewidth = 0.1,alpha = 0.4)  + 
  geom_point(size=2.5, aes(color =Dosis ,  shape=Treatment),stroke = 1) + #size = 4,
  labs(x='PCoA1', y='PCoA2')+
  theme_few()
ylim(-0.1,0.15)

f
ggplot(data = PCoA_Comp_mean2,aes(mean.x,mean.y)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),linewidth = 0.1,alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),linewidth = 0.1,alpha = 0.4)  + 
  geom_point(size=2.5, aes(color =Dosis ,  shape=Treatment),stroke = 1) + #size = 4,
  labs(x='PCoA1', y='PCoA2')+
  theme_few()
ylim(-0.1,0.15)

