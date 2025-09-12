###################################################################
library("FactoMineR")
library("factoextra")
library(ade4)
library(MASS)
library(ellipse)
library(ggplot2)
library(agricolae)
library(vegan)
library(ggthemes)
library(RVAideMemoire)
library(reshape2)
library(emmeans)
library(dplyr)
library(pals)
library(paletteer)

setwd("G:/My Drive/labs/Nottingham/Duckweed/exp Metabolites/Metabolome/")

Metabolome <- read.table("metabolite_log2foldchange2", header = TRUE, row.names = 1)

Species=c('SP7498','SI7820','LT9243','LA7339','SP9509','LY9205','SI9227','LJ9250','LP0049')
Complexity=c(7.25,6.8,6.2,5,7,6.8,7.6,5.33,5)
Met=data.frame(cbind(Species,Complexity))

TableMet6 = merge(list_log2foldchange2, Met, by= 'Species')

unique(mms_all2$Metabolite) 

mms_all2$log2foldchange = as.numeric(mms_all2$log2foldchange)
mms_all2$Complexity = as.numeric(mms_all2$Complexity)

paleta_alive <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                  "#661100","#44AA99", "#999933", "#882255",  "#6699CC", "#888888",
                  "#49AA91", "#999983", "#882249", "#661830", "#6699DC", "#828389")


corrX2262 = ggplot(TableMet6[TableMet6$Metabolite=='X2322',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method=lm, se=TRUE)+
  theme(legend.position='none')+
  labs(y='log2foldchange', x='CCL',
       title='X2322')+
  scale_color_manual(values = paleta_alive)

summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X2262',]))
#Multiple R-squared:  0.4508,	Adjusted R-squared:  0.3723 
#F-statistic: 5.746 on 1 and 7 DF,  p-value: 0.04767


corrX1942 = ggplot(mms_all2[mms_all2$Metabolite=='X1942',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method=lm, se=TRUE)+
  theme(legend.position='none')+
  labs(y='log2foldchange', x='CCL',
       title='X1942')+
  scale_color_manual(values = paleta_alive)


summary(lm(Complexity~log2foldchange, data=mms_all2[mms_all2$Metabolite=='X1942',]))
#Multiple R-squared:  0.4933,	Adjusted R-squared:  0.4209 
#F-statistic: 6.815 on 1 and 7 DF,  p-value: 0.03489


corrX2175 = ggplot(mms_all2[mms_all2$Metabolite=='X2175',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position ='none' )+
  labs(y='log2foldchange', x='CCL',
       title='X2175')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X2175',]))
#Multiple R-squared:  0.4504,	Adjusted R-squared:  0.3719 
#F-statistic: 5.736 on 1 and 7 DF,  p-value: 0.04781


corrX1624 = ggplot(mms_all2[mms_all2$Metabolite=='X1624',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X1624')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X1624',]))
#Multiple R-squared:  0.6091,	Adjusted R-squared:  0.5533 
#F-statistic: 10.91 on 1 and 7 DF,  p-value: 0.01307


corrX2290 = ggplot(mms_all2[mms_all2$Metabolite=='X2290',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X2290')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X2290',]))
#Multiple R-squared:  0.6635,	Adjusted R-squared:  0.6155  
#F-statistic: 13.8  on 1 and 7 DF,  p-value: 0.0075



corrX679 = ggplot(mms_all2[mms_all2$Metabolite=='X679',], #
                  aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X679')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X679',]))
#Multiple R-squared:  0.5688,	Adjusted R-squared:  0.5072 
#F-statistic: 9.235 on 1 and 7 DF,  p-value: 0.01887



corrX1290 = ggplot(mms_all2[mms_all2$Metabolite=='X1290',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X1290')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X1290',]))
#Multiple R-squared:  0.5273,	Adjusted R-squared:  0.4598 
#F-statistic: 7.808 on 1 and 7 DF,  p-value: 0.02674

corrX1221 = ggplot(mms_all2[mms_all2$Metabolite=='X1221',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X1221')+
  scale_color_manual(values = paleta_alive)

summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X1221',]))
#Multiple R-squared:  0.4876,	Adjusted R-squared:  0.4144 
#F-statistic:  6.66 on 1 and 7 DF,  p-value: 0.03643


corrX712 = ggplot(mms_all2[mms_all2$Metabolite=='X712',], #
                  aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position='none')+
  labs(y='log2foldchange', x='CCL',
       title='X712')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X712',]))
#Multiple R-squared:  0.5352,	Adjusted R-squared:  0.4688 
#F-statistic: 8.059 on 1 and 7 DF,  p-value: 0.02508


corrX1082 = ggplot(mms_all2[mms_all2$Metabolite=='X1082',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X1082')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X1082',]))
#Multiple R-squared:  0.4467,	Adjusted R-squared:  0.3676 
#F-statistic:  5.65 on 1 and 7 DF,  p-value: 0.04909

corrX182 = ggplot(mms_all2[mms_all2$Metabolite=='X182',], #
                  aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X182')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X182',]))
#Multiple R-squared:  0.5112,	Adjusted R-squared:  0.4414 
#F-statistic: 7.322 on 1 and 7 DF,  p-value: 0.03038


corrX1454 = ggplot(mms_all2[mms_all2$Metabolite=='X1454',], #
                   aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X1454')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X1454',]))
#Multiple R-squared:  0.7419,	Adjusted R-squared:  0.7051 
#F-statistic: 20.13 on 1 and 7 DF,  p-value: 0.002845


corrX358 = ggplot(mms_all2[mms_all2$Metabolite=='X358',], #
                  aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='358')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X358',]))
#Multiple R-squared:  0.6434,	Adjusted R-squared:  0.5925 
#F-statistic: 12.63 on 1 and 7 DF,  p-value: 0.009294


corrX984 = ggplot(mms_all2[mms_all2$Metabolite=='X984',], #
                  aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X984')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X984',]))
#Multiple R-squared:  0.4994,	Adjusted R-squared:  0.4279 
#F-statistic: 6.983 on 1 and 7 DF,  p-value: 0.0333


corrX822 = ggplot(mms_all2[mms_all2$Metabolite=='X822',], #
                  aes(x=Complexity, y=log2foldchange))+
  theme_bw() + 
  geom_point(size=3, aes(col=Species)) + 
  geom_smooth(method="lm", se=T)+
  theme(legend.position = 'none')+
  labs(y='log2foldchange', x='CCL',
       title='X822')+
  scale_color_manual(values = paleta_alive)


summary(lm(log2foldchange~Complexity, data=mms_all2[mms_all2$Metabolite=='X822',]))
#Multiple R-squared:  0.541,	Adjusted R-squared:  0.4754 
#F-statistic:  8.25 on 1 and 7 DF,  p-value: 0.02391

library(ggpubr)
ggarrange(corrX1942, corrX2175, corrX1624, corrX2290,
          corrX679, corrX1290, corrX712, corrX2262,
          corrX1082, corrX182, corrX1454, corrX358,
          corrX984, corrX822, corrX1221,  
          ncol = 4, nrow=4)
