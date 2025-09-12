###Heatmap
library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(ggdendro)
library(factoextra)
library(dplyr)
library(ggdendro)
library(reshape2)
library(cowplot)

setwd("G:/My Drive/labs/Nottingham/Duckweed/thermodynamics/")

###Data
Data <- read.table("Themodynamics.txt", header = TRUE, row.names = 1)
Data[2,1:15]

Data.mean <- Data %>% group_by(Species, Treatment,Root_complexity )%>%
  summarise_all(mean)


# Scatterplot
Data2 = Data[-c(1,54),]
Data2$Root_complexity = as.numeric(Data2$Root_complexity)
theme_set(theme_bw())  # pre-set the bw theme.

ggplot(Data2[Data2$Treatment == 'NB',], aes(Root_complexity, Fmax)) +  
  geom_jitter(aes(col=Species), size = 2.5) + 
  geom_smooth(method="lm", se=T, col = 'black')+
  theme_few()+
  labs(y='Maximun peak', x= 'CCL',col='black')
