####Experiemnt2
library(ggplot2)
library(ggthemes)
library(ggpubr)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/")
Growth <- read.table("Exp2_Duckweed_Curve_Growth.txt", header = TRUE, row.names = 1)
Growth$Nutrient = factor(Growth$Nutrient, c('High','Medium','Low'))

SI7820 = ggplot(Growth[Growth$Plant == '7820',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment))  +   # , pch=Nutrient draw points
  facet_grid(.~Nutrient, scales = "fixed",space = "fixed")+
  geom_smooth(aes(col=Treatment, linetype=Nutrient),method="loess", se=T, level = 0.8) +
  labs(title='SI7820',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="none")

SP9509 = ggplot(Growth[Growth$Plant == '9509',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment))  +   #, pch=Nutrient draw points
  facet_grid(.~Nutrient, scales = "fixed",space = "fixed")+
  geom_smooth(aes(col=Treatment),method="loess", se=T, level = 0.8) +
  labs(title='SP9509',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="  none")+
  ylim(1,3)

LP8539 = ggplot(Growth[Growth$Plant == '8539',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment))  +   # draw points
  facet_grid(.~Nutrient, scales = "fixed",space = "fixed")+
  geom_smooth(aes(col=Treatment),method="loess", se=T, level = 0.8) +
  labs(title='LP8539',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="  none")+
  ylim(1,3)

LA7339 = ggplot(Growth[Growth$Plant == '7339',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment))  +   # draw points
  facet_grid(.~Nutrient, scales = "fixed",space = "fixed")+
  geom_smooth(aes(col=Treatment),method="loess", se=T, level = 0.8) +
  labs(title='LA7339',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="  none")

PS = ggplot(Growth[Growth$Plant == 'PS',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment, pch=Nutrient))  +   # draw points
  facet_grid(.~Nutrient, scales = "fixed",space = "fixed")+
  geom_smooth(aes(col=Treatment),method="loess", se=T, level = 0.8) +
  labs(title='PS',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="  none")
#LP0049
LP0049 = ggplot(Growth[Growth$Plant == '0049',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment))  +   # draw points
  facet_grid(.~Nutrient, scales = "fixed",space = "fixed")+
  geom_smooth(aes(col=Treatment),method="loess", se=T, level = 0.8) +
  labs(title='LP0049',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="  none")+
  ylim(1,4)

#LJ92509250
LJ9250 = ggplot(Growth[Growth$Plant == '9250',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment))  +   # draw points
  facet_grid(.~Nutrient, scales = "fixed",space = "fixed")+
  geom_smooth(aes(col=Treatment),method="loess", se=T, level = 0.8) +
  labs(title='LJ9250',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="  none")

ggplot(Growth[Growth$Plant == '9250',], aes(x=Day, y=Value, fill=MxT)) + 
  geom_point(aes(col=Treatment))  +   # draw points
  geom_smooth(aes(col=Treatment, linetype=Nutrient),method="loess", se=F, level = 0.8) +
  labs(title='LJ9250',  y='cm2', size = 0.4) +
  theme_few()+ theme(legend.position="top")

ggarrange(SP9509, SI7820,LP0049,LJ9250,LA7339,LP8539, 
          ncol = 1, nrow=6)#PS,
