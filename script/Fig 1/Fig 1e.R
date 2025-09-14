library(ggplot2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp1/")

Status <- read.table("Status.txt", header = TRUE, row.names = 1)

Status$Species <- factor(Status$Species,rev(c('PS','SP7820','SI9227','SP9509','SP5543',
                                              'SP7498','SP9192','LP7760','LP0049','LT9109',
                                              'LJ9250','LM7123','LM8389','LV7650','LY5280',
                                              'LM7200','LA7339','LP8539','LT9243')))
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
 
