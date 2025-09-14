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

setwd("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 4/Fig 4a/")
Exudates <- read.table("Exudates.txt", header = TRUE, row.names = 1)

Exudates[,1:5]

varespec.bray <- vegdist(Exudates[,-c(1:4)], method = "bray") # -c(3,27)dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan
pcoaVS <- pcoa(varespec.bray)
pcoaVS$vectors[,1:2]

f<-data.frame(pcoaVS$vectors[,1:4])
f <- cbind(Exudates[,c(1:4)],f )#-c(3,27)
rownames(f)<-x$Row.names
f$Row.names <- NULL

PCA_Comp_mean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2,CCL=CCL,CCL_ds=CCL_ds)~Species ,f,mean)
PCA_Comp_meands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Species,f,sd)
PCoA_Comp_mean2 =cbind(PCA_Comp_mean,PCA_Comp_meands[,c(2:3)])


paleta_alive <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                  "#661100","#44AA99", "#999933", "#882255",  "#6699CC", "#888888",
                  "#49AA91", "#999983", "#882249", "#661830", "#6699DC", "#828389")

ggplot(data = PCoA_Comp_mean2, aes(mean.x,mean.y)) + #PCoA_Comp.F_mean2 mean.x mean.y
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),linewidth = 0.1,alpha = 0.3) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),linewidth = 0.01,alpha = 0.1)  + 
  geom_point(size=2.5, aes(color = Species),stroke = 1) + #size = 4,
  labs(x='PCoA1', y='PCoA2') +
  #scale_color_gradientn(colours=c('blue','red'))+
  theme_few()+
  scale_color_manual(values = paleta_alive)

#####PERMANOVA
###PERMANOVA
# Calculate bray curtis distance matrix
ps.3_bray <- phyloseq::distance(ps.Syncom.Root5, method = "bray")
# make a data frame from the sample_data
sample_ps.3c <- data.frame(sample_data(ps.Syncom.Root5))
# Adonis test
adonis2(varespec.bray ~ CCL*Species, data = Exudates)


