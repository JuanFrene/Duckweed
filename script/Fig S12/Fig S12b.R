library(ggplot2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

Status <- read.table("Plant fitness.txt", header = TRUE, row.names = 1)
Status$Microbiome[Status$Microbiome=="NB"]<-"Sterile Medium"
Status$Compound[Status$Compound=="Glutaric_acid"]<-"Glutaric acid"

Status.mean <- Status %>% group_by(Species,Microbiome, Compound)%>%
  summarise_all(mean)

melted_Status <- Status.mean[,-c(4:9)] %>% melt
melted_Status$variable <- factor(melted_Status$variable,c('Dead','Stressed','Healthy'))
melted_Status$Compound <- factor(melted_Status$Compound,c('None','DMSO','Methyl-Cytrate','Citrulline','Adrenaline','Lysine','Glutaric acid','Pyrone'))

melted_Status2 = melted_Status[melted_Status$Microbiome != 'Sterile Medium',]
### Fitness figure Fitness bar plot
ggplot(data=melted_Status2[melted_Status2$Compound != 'DMSO',], aes(y=Species, x=value, fill=variable)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  theme(legend.position="top") +
  geom_bar(stat="identity", color = "black") + 
  facet_grid(~Compound, space = "free",scales = "free") +
  theme_few() +
  scale_fill_manual(values=c("#C1CDCD", "#FFD700","#66CD00"))+ 
  labs(y='Species', x='Porcentage (%)') 
