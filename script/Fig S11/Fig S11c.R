library(ggplot2)
library(paletteer)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp1/")

Model <- read.table( "Growth and vel.txt", header = TRUE, row.names = 1)
Root_Trait <- read.table( "Root traits.txt", header = TRUE, row.names = 1)
head(Model)
Model$Species <- factor(Model$Species,rev(c('PS','SP7820','SI9227','SP9509','SP5543',
                                      'SP7498','SP9192','LP7760','LP0049','LT9109',
                                      'LJ9250','LM7123','LM8389','LV7650','LY5280',
                                      'LM7200','LA7339','LP8539','LT9243')))



#############
correlation
Model ##Growth parameter table
Root_Trait ##Root trait parameter table

nrow(Root_Trait)
nrow(Model)
All = cbind(Root_Trait, Model[,-(2:4)])
ncol(Root_Trait)

##Corr y p-value
CorrNB<-cor(All[1:95,-c(1:3,25:35)], method = 'pearson')
p.matnB <-cor_pmat(All[1:95,-c(1:3,25:35)], method = 'pearson')

CorrB<-cor(All[-(1:95),-c(1:3,25:35)], method = 'pearson')
p.matB <-cor_pmat(All[-(1:95),-c(1:3,25:35)], method = 'pearson')


in_vel = data.frame(cbind(rownames(CorrB), CorrB[,24], p.matB[,24]))
colnames(in_vel) = c('Root','r', 'p.value')
in_vel$r = as.numeric(in_vel$r)

####Adjust the p-value
in_vel$SignPvalue <- "NoSignificant"
pval_thres <- 0.05
in_vel$SignPvalue[which(in_vel$p.value < pval_thres)] <- "p < 0.05"
in_vel$SignPvalue <- in_vel$SignPvalue %>% factor

ggplot(in_vel[-(22:27),], aes(r,Root)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_few()+ 
  theme(legend.position="top")+
  geom_tile(aes(color = SignPvalue),fill = '#00000000', linewidth = 0.2,width = 0.01,height = 0.95) +
  xlab('r')+
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significative")
  