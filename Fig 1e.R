library(ggplot2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

Status <- read.table("Plant fitness.txt", header = TRUE, row.names = 1)
Status$Species <- factor(Status$Species,rev(c('LY5280','SP7498','SP9509',
                                              'LJ9250','LP0049','LP8539')))
Status$Fitness <- factor(Status$Fitness,c('Dead','Stressed','Healthy'))
Status$Treatment[Status$Treatment=="No-Microbiome"]<-"Sterile Water"
Status$Treatment[Status$Treatment=="With-Microbiome"]<-"Natural Water"

Status$Treatment <- factor(Status$Treatment,c('Sterile Water','Natural Water'))

### Fitness figure
ggplot(data=Status, aes(y=Species, x=Percentage, fill=Fitness)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  theme(legend.position="top") +
  geom_bar(stat="identity", color = "black") + 
  facet_grid(~Treatment, space = "free",scales = "free") +
  theme_few() +
  scale_fill_manual(values=c("#C1CDCD", "#FFD700","#66CD00"))+ 
  labs( y='Species') 
 
