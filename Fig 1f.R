library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(harrietr)
library(egg)
library(paletteer)
library(scales)
library(car)
library(Rmisc)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/Exp1/")

TablaCompleta = read.table("Root traits.txt", header = TRUE, row.names = 1)

Mean = aggregate(cbind(mean.x=cell_layer_n,mean.y=Dead)~Species*Microbiome,TablaCompleta,mean)
SD = aggregate(cbind(ds.x=cell_layer_n,ds.y=Dead)~Species*Microbiome,TablaCompleta,se)
Table2 <- cbind(Mean,SD[3:4])

ggplot(data = Table2, aes(mean.x, mean.y, )) +
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.2, alpha = 0.4) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),size = 0.2, alpha = 0.4)  + 
  geom_hline(yintercept = 0,size = 0.2,color = "black",linetype = "longdash")+
  geom_point( size = 2, aes(colour = Microbiome),stroke = 1) + 
  geom_smooth(method="lm",size = 0.5, se=F,aes(linetype = Microbiome),color='black') + 
  labs(y='Mortality (%)', x='CCL') +
  theme_few()+
  ylim(-20,100)+
  scale_colour_manual(values=c('orange', '#0000ff'),l=c('Natural Water','Sterile Water'))
  
summary(lm(cell_layer_n~Dead, data = TablaCompleta[TablaCompleta$Microbiome=='NB',]))
#NB Multiple R-squared:  0.004726 p-value: 0.508

summary(lm(cell_layer_n~Dead, data = TablaCompleta[TablaCompleta$Microbiome!='NB',]))
#NB Multiple R-squared:  0.1576 p-value: 6.804e-05


