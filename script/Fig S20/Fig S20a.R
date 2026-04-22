######Load packages
library(Rmisc)
library(reshape2)
library(ggrepel)
library(scales)
library(ggtree)
library(harrietr)
library(emmeans)
library(paletteer)
library(pals)
library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(ggdendro)



######Set work directory
setwd("G:/My Drive/labs/Nottingham/Duckweed/Drop out/Analysis files/")

######Open data file
Data <- read.table("SP9509 dropout.txt", header = TRUE, row.names = 1)
head(Data)
ID=rownames(Data)
Data= cbind(ID,Data)
Data$Syncom[which(Data$Syncom == 'Full_Syncom')] <- "Syncom"

#Barplot complexity differences
display_D2$complexity
display_D3 = cbind(rownames(display_D2),display_D2)
colnames(display_D3) = c('Syncom', colnames(display_D2))
Table = Res_p_Drop[Res_p_Drop$UId == '_Complexity',]
Table$Syncom = factor(Table$Syncom, c('R245-Syncom','R289-Syncom',"NB-Syncom","P9-Syncom",'P153-Syncom','P234-Syncom',"P162-Syncom","P129-Syncom","P118A-Syncom","P128A.22-Syncom","R34-Syncom","R148A-Syncom","R345.2-Syncom","R186-Syncom","P74A-Syncom","R151-Syncom","R28-Syncom","P69-Syncom","P192-Syncom","R88-Syncom","P65-Syncom"))
head(Res_p_Drop[Res_p_Drop$UId == '_Complexity',])

ggplot(data=Table[-17,], aes(y=Syncom, x=t.ratio)) +
  geom_bar(stat="identity", colour='black')+
  theme_few()+ labs(x = 'Root traits',y = 'Syncom')+
  coord_flip()+ theme(legend.position = 'none')


