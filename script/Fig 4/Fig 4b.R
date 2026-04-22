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

setwd("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 4/Fig 4b/")
Metabolome <- read.table("Table_metabolome.txt", header = TRUE, row.names = 1)
Metabolome1 = Metabolome[Metabolome$Treatment != 'QC',]
Metabolome2 = Metabolome1[Metabolome1$Species != 'LP8539',]

Metabolome1[,1:7]

varespec.bray <- vegdist(Metabolome2[,-c(1:5)], method = "euclidea") # -c(3,27)dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan
pcoaVS <- pcoa(varespec.bray)
pcoaVS$vectors[,1:2]

f<-data.frame(pcoaVS$vectors[,1:2])
f <- cbind(Metabolome2[,c(1:5)],f )#-c(3,27)
rownames(f)<-f$Row.names
f$Row.names <- NULL

PCA_Comp_mean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2, CCL=Complexity)~Species*Treatment ,f,mean)
PCA_Comp_meands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Species*Treatment,f,sd)
PCoA_Comp_mean2 =cbind(PCA_Comp_mean,PCA_Comp_meands[,c(3:4)])
PCoA_Comp_mean2$Treatment
ggplot(data = Table, aes(CCL, PCoA1)) +
  #geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),linewidth = 0.2, alpha = 0.4) +
  #geom_errorbar(mapping = aes(ymin =mean.y - ds.x ,ymax = mean.x + ds.x),linewidth = 0.2, alpha = 0.4)  + 
  geom_hline(yintercept = 0,size = 0.2,color = "black",linetype = "longdash")+
  facet_grid(.~Treatment)+
  geom_point( size = 2, aes(colour = Species),stroke = 1) + 
  geom_smooth(method="lm",size = 0.5, se=T,,color='black') + 
  labs(y='PCoA 1', x='CCL') +
  theme_few()
  scale_color_manual(values = paleta_alive)
  
summary(lm(CCL~PCoA1, data=Table[Table$Treatment=='NB',]))  
#Multiple R-squared:  0.4035,	Adjusted R-squared:  0.3182 
#F-statistic: 4.734 on 1 and 7 DF,  p-value: 0.06605

summary(lm(PCoA1~CCL, data=Table[Table$Treatment!='NB',]))  
#Multiple R-squared:  0.4474,	Adjusted R-squared:  0.3684 
#F-statistic: 5.667 on 1 and 7 DF,  p-value: 0.04884

Species= c("LA7339","LJ9250","LP0049", "LT9243", "LY9205", "SI7820", "SI9227",
        "SP7498", "SP9509", "LA7339", "LJ9250", "LP0049" ,"LT9243" ,"LY9205",
        "SI7820", "SI9227","SP7498","SP9509")
Treatment=c("NB","NB","NB","NB","NB","NB","NB","NB","NB","Syncom","Syncom", "Syncom", "Syncom", "Syncom",
          "Syncom","Syncom","Syncom","Syncom")

CCL = c(5,5,5,5.3938461538461535, 5.793846153846154,5.787692307692307,7.313846153846153,6.803076923076922,6.587692307692308,4.993846153846153,
        5.332307692307692,5.196923076923077,6.193846153846153,7.603076923076923,6.796923076923077,7.596923076923076,7.252307692307692,6.993846153846153)
        

PCoA1 = c(15.480106100795759,34.23010610079576, 9.661140583554378, 8.32559681697613, 1.6014588859416463, -5.294429708222808, -2.4416445623342256, 
          -0.44694960212201806, 7.765915119363385, -1.329575596816973, 25.573607427055702, 5.976127320954909, -9.64854111405836, -11.524535809018573,
          -27.816976127320956, -9.368700265251995, -29.590185676392565, -6.717506631299742)

Table =data.frame(cbind(Species, Treatment,CCL,PCoA1))
Table$CCL =as.numeric(Table$CCL)
Table$PCoA1 =as.numeric(Table$PCoA1)
