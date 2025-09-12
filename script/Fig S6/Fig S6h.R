library(devtools)
install_github("tomhopper/windrose")
library(windrose)
library(agricolae)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/")

TablaCompleta = read.table("comparison exp1 exp2.txt", header = TRUE, row.names = 1)

TablaCompleta$Species = factor(TablaCompleta$Species, c('SP9509','SI7820','LJ9250','LP0049','LA7339','LP8539','LY5280')) 

TablaCompletaL = TablaCompleta[TablaCompleta$Nutrition=='Low',]
TablaCompletaM = TablaCompleta[TablaCompleta$Species=='Medium',]
TablaCompletaH = TablaCompleta[TablaCompleta$Experiment=='High',]

########## Bar plot y Polar bar plot de los root traits 
#Tables
Coordata <- read.table("clipboard", header = TRUE)
str(Coordata)
Coordata$Microbiome = factor(Coordata$Microbiome, c('SynCom', 'NB'))
Coordata$Cell_type = factor(Coordata$Cell_type, c('Epidermis','Exodermis','CCL4','CCL3','CCL2','CCL1','Endodermis','Stele')) 
Coordata$Species = factor(Coordata$Species, c('LY5280','SP9509','SI7820','LP0049','LJ9250','LA7339','LP8539')) 

Coordata2=Coordata[Coordata$Experiment=='1',]
CoordataV2 = Coordata[Coordata$Nutrition!='1.34mM',]
CoordataV3 = CoordataV2[CoordataV2$Species!='LY5280',]

Mean.Traits <- CoordataV3%>%group_by(Nutrition,Microbiome,Species,Cell_type, )%>%
  summarise_all(mean)

Mean.TraitsB = Mean.Traits[Mean.Traits$Microbiome!='NB',]
Mean.TraitsBHigh = Mean.TraitsB[Mean.TraitsB$Nutrition=='335mM',]
Mean.TraitsBMedium = Mean.TraitsB[Mean.TraitsB$Nutrition=='6.7mM',]
Mean.TraitsBLow = Mean.TraitsB[Mean.TraitsB$Nutrition=='0.67mM',]

Mean.TraitsNB = Mean.Traits[Mean.Traits$Microbiome=='NB',]
Mean.TraitsNBHigh = Mean.TraitsNB[Mean.TraitsNB$Nutrition=='335mM',]
Mean.TraitsNBMedium = Mean.TraitsNB[Mean.TraitsNB$Nutrition=='6.7mM',]
Mean.TraitsNBLow = Mean.TraitsNB[Mean.TraitsNB$Nutrition=='0.67mM',]

pBHigh = ggplot(data=Mean.TraitsBHigh, aes(x=Species, y=Diameter, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #theme_few()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") + 
  ylim(0,100)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pBCircHigh = pBHigh + coord_polar()+ theme(legend.position="none")+xlab('')

pBMedium = ggplot(data=Mean.TraitsBMedium, aes(x=Species, y=Diameter, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #theme_few()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,100)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pBCircMedium = pBMedium + coord_polar()+ theme(legend.position="none")+xlab('')

pBLow = ggplot(data=Mean.TraitsBLow, aes(x=Species, y=Diameter, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #theme_few()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,100)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pBCircLow = pBLow + coord_polar()+ theme(legend.position="none")+xlab('')

pNBHigh = ggplot(data=Mean.TraitsNBHigh, aes(x=Species, y=Diameter, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #theme_few()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,100)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pNBCircHigh = pNBHigh + coord_polar()+ theme(legend.position="none")+xlab('')

pNBMedium = ggplot(data=Mean.TraitsNBMedium, aes(x=Species, y=Diameter, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #theme_few()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,100)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pNBCircMedium = pNBMedium + coord_polar()+ theme(legend.position="none")+xlab('')

pNBLow = ggplot(data=Mean.TraitsNBLow, aes(x=Species, y=Diameter, fill=Cell_type)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #theme_few()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") + ylim(0,100) + 
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pNBCircLow = pNBLow + coord_polar()+ theme(legend.position="none")+xlab('')

ggarrange(pBCircHigh,pNBCircHigh,pBCircMedium,pNBCircMedium,pBCircLow,pNBCircLow, ncol = 2, nrow=3)



g = ggplot(data=Mean.TraitsB, aes(x=Species, y=Porcentual, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+theme_few()+ theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
gBCirc = g + coord_polar() + theme(axis.text.x = element_text(angle = 45, hjust=1),axis.text.y = element_blank())+ylab('')+xlab('')


g = ggplot(data=Mean.TraitsNB, aes(x=Species, y=Porcentual, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+theme_few()+ theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
gNBCirc = g + coord_polar() + theme(axis.text.x = element_text(angle = 45, hjust=1),axis.text.y = element_blank())+ylab('')+xlab('')

ggarrange(gNBCirc,gBCirc,
          ncol = 2, nrow=1)

