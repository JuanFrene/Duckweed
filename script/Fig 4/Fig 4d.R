library(ggplot2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

Status <- read.table("Plant fitness LC50.txt", header = TRUE, row.names = 1)
Status$Microbiome[Status$Microbiome=="NB"]<-"Sterile Medium"

Status.mean <- Status %>% group_by(Microbiome,Compound,Dosis,ID)%>%
  summarise_all(mean)

melted_Status <- Status.mean[,-c(5:9)] %>% melt
melted_Status$variable <- factor(melted_Status$variable,c('Dead','Stressed','Healthy'))

melted_Status$Compound <- factor(melted_Status$Compound,c("NB","Syncom","C4","C5"))
melted_Status$Dosis <- factor(melted_Status$Dosis,rev(c("NB","Syncom","25", "50","75")))

### Fitness figure Fitness bar plot
ggplot(data=melted_Status[melted_Status$ID!='Syncom',], aes(y=Dosis, x=value, fill=variable)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  theme(legend.position="top") +
  geom_bar(stat="identity", color = "black") + 
  facet_grid(Compound~., space = "free",scales = "free") +
  theme_few() +
  scale_fill_manual(values=c("#C1CDCD", "#FFD700","#66CD00"))+ 
  labs(y='Species', x='Porcentage (%)') 
 
