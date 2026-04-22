library(ggplot2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/")

Status <- read.table("Exp2_Duckweed_fitness.txt", header = TRUE, row.names = 1)

Status$Species <- factor(Status$Species,rev(c('SP9509','SI7820','LJ9250','LP0049','LA7339','LP8539')))
Status$Fitness <- factor(Status$Fitness,c('Dead','Stressed','Healthy'))
Status$Treatment[Status$Treatment=="NB"]<-"Sterile Water"

Status$Treatment <- factor(Status$Treatment,c('Sterile Water','Syncom'))

Fitness = Status%>%group_by(Species, Treatment, Fitness, Medium)%>%
  summarise_all(mean)
Fitness = data.frame(Fitness)

### Fitness figure
ggplot(data=Fitness, aes(y=Species, x=Value, fill=Fitness)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  theme(legend.position="top") +
  geom_bar(stat="identity", color = "black") + 
  facet_grid(Medium~Treatment, space = "free",scales = "free") +
  theme_few() +
  scale_fill_manual(values=c("#C1CDCD", "#FFD700","#66CD00"))+ 
  labs( y='Species', x = 'Porcentage') 

