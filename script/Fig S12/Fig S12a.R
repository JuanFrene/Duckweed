packages <- c('ggcorrplot','ggthemes','dplyr', "ape",  "Biostrings","tibble", 
              "lme4", "lmerTest", "ggplot2", "vegan", "car", "rcompanion",'ggpubr', 
              "emmeans", "RVAideMemoire",'gplots','plotly','tidyr', 'reshape2')
sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

Growth2 <- read.table("Plant Growth.txt", header = TRUE, row.names = 1)
colnames(Growth2) = c('Species', 'Microbiome','Compound', 'Placa', 'Well', '0', '3','6','12','15','18','21')

Growth2$Placa = as.factor(Growth2$Placa)
Growth2$Well = as.factor(Growth2$Well)

melt.Species <- melt(Growth2)
melt.Species$value = as.numeric(melt.Species$value)
melt.Species$variable = as.numeric(melt.Species$variable)
melt.Species2 <- melt.Species[melt.Species$Microbiome=='Syncom',]
melt.Species3 <- melt.Species2[melt.Species2$Compound!='DMSO',]

#Growth curve
paleta_alive <- c("#C200FF",'#FF0000','#8B7500','#00008B',"#FFB919",'#FF7F50',"#00CC7A", 'black', 'grey')

LY5290 <- ggplot(melt.Species3[melt.Species3$Species=='LY5290',], aes(x=variable, y=value, fill=Compound)) + 
  geom_point(aes(col=Compound)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8, aes(col=Compound)) + 
  labs(title='LY5290',  y='cm2', x= 'Days')+
  scale_color_manual(values = paleta_alive)

LJ9250 <- ggplot(melt.Species3[melt.Species3$Species=='LJ9250',], aes(x=variable, y=value, fill=Compound)) + 
  geom_point(aes(col=Compound)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8, aes(col=Compound)) + 
  labs(title='LJ9250',  y='cm2', x= 'Days')+
  scale_color_manual(values = paleta_alive)

SP9509 <- ggplot(melt.Species3[melt.Species3$Species=='SP9509',], aes(x=variable, y=value, fill=Compound)) + 
  geom_point(aes(col=Compound)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8, aes(col=Compound)) + 
  labs(title='SP9509',  y='cm2', x= 'Days')+
  scale_color_manual(values = paleta_alive)

LP8539 <- ggplot(melt.Species3[melt.Species3$Species=='LP8539',], aes(x=variable, y=value, fill=Compound)) + 
  geom_point(aes(col=Compound)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8, aes(col=Compound)) + 
  labs(title='LP8539',  y='cm2', x= 'Days')+
  scale_color_manual(values = paleta_alive)

LP0049 <- ggplot(melt.Species3[melt.Species3$Species=='LP0049',], aes(x=variable, y=value, fill=Compound)) + 
  geom_point(aes(col=Compound)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8, aes(col=Compound)) + 
  labs(title='LP0049',  y='cm2', x= 'Days')+
  scale_color_manual(values = paleta_alive)

SP7498 <- ggplot(melt.Species3[melt.Species3$Species=='SP7498',], aes(x=variable, y=value, fill=Compound)) + 
  geom_point(aes(col=Compound)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8, aes(col=Compound)) + 
  labs(title='SP7498',  y='cm2', x= 'Days')+
  scale_color_manual(values = paleta_alive)

ggarrange(LY5290,SP7498,SP9509,
          LJ9250,LP0049,LP8539,common.legend = TRUE, legend = 'right', 
          ncol = 2, nrow=3)



