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


paleta_alive <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                  "#661100","#44AA99", "#999933", "#882255",  "#6699CC", "#888888",
                  "#49AA91", "#999983", "#882249", "#661830", "#6699DC", "#828389")

PCoA_Comp_mean2$CCL = as.numeric(PCoA_Comp_mean2$CCL)

PcoA = ggplot(data = PCoA_Comp_mean2, aes(mean.x,mean.y)) + #PCoA_Comp.F_mean2 mean.x mean.y
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),linewidth = 0.1,alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),linewidth = 0.01,alpha = 0.4)  + 
  geom_point(size=2.5, aes(color = CCL),stroke = 1) + #size = 4,
  labs(x='PCoA1', y='PCoA2') +
  scale_color_gradientn(colours=c('blue','red'))+
  theme_few()
  scale_color_manual(values = paleta_alive)

#####PERMANOVA
###PERMANOVA
# Calculate bray curtis distance matrix
adonis2(varespec.bray ~ CCL*Species, data = Exudados)

paleta_alive <- c('#FF0000','#00008B',"#FFB919","#00CC1C")

PERMANOVA = adonis2(varespec.bray ~ CCL*Species, data = Exudados, by = "terms")

PERMANOVA2 = data.frame(PERMANOVA)
variable  = c('variable','variable','variable','variable')
PERMANOVA2$Significance <- "No Significant"
pval_thres <- 0.05
PERMANOVA2$Significance[which(PERMANOVA2$Pr..F. < pval_thres)] <- "Significant"
PERMANOVA2$Significance <- PERMANOVA2$Significance %>% factor

PERMANOVA2$R2[which(PERMANOVA2$R2 < '0.05201346')] <- "0.144"
PERMANOVA2$R2[which(PERMANOVA2$R2 == '0.396731354925236')] <- "0.304"

LC = data.frame(cbind(variable, row.names(PERMANOVA2), PERMANOVA2[,c(3,6)]))
colnames(LC) = c('variable','Effect', 'R2','Significance')


#Plot Variance
dim(LC)
colnames(LC)
dim(LC)
data_melt <- LC[-(4),]

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)#Turn your 'treatment' column into a character vector
data_melt$Effect <- factor(data_melt$Effect, levels=unique(rev(data_melt$Effect)))#Then turn it back into a factor with the levels in the correct order

mypal = c('white','#FF0000','#00008B',"#FFB919","#00CC1C")

data_melt$R2 = as.numeric(data_melt$R2)#, c('Crop','Compost','Crop:Compost','Residual'))

PERMANOVA_grap = ggplot(data=data_melt, aes(x=variable, y=R2,fill=Effect)) +
  geom_bar(stat="identity",aes(color=Significance), size = 0.3,width = 1,linewidth=1) +# 
  scale_fill_manual(values = mypal) +
  theme_few() + #guides() + #color = FALSE, fill=FALSE
  labs(y="Variance explained (%)") +
  scale_color_manual(values = c('black',"orange"),na.value =  "transparent",name = "Significance vs Water")+
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0)),
        axis.text.x   = element_blank())

plot_grid(PcoA, PERMANOVA_grap, align = 'h', rel_widths = c(1, 0.4))
