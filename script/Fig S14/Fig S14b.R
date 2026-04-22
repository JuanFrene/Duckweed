library(devtools)
#install_github("tomhopper/windrose")
library(windrose)
library(agricolae)
library(reshape2)
library(dplyr)
packages <- c('ggthemes', "ape", "ShortRead", "Biostrings", "phyloseq", 
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4", 
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire")
sapply(packages, require, character.only = TRUE)              


setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/")

TablaCompleta = read.table("comparison exp1 exp2.txt", header = TRUE, row.names = 1)

TablaCompleta$Species = factor(TablaCompleta$Species, c('SP9509','SI7820','LJ9250','LP0049','LA7339','LP8539','LY5280')) 

TablaCompletaL = TablaCompleta[TablaCompleta$Nutrition=='Low',]
TablaCompletaM = TablaCompleta[TablaCompleta$Species=='Medium',]
TablaCompletaH = TablaCompleta[TablaCompleta$Experiment=='High',]

########## Bar plot y Polar bar plot de los root traits 
#Tables
Coordata <- read.table("G:/My Drive/labs/Nottingham/Duckweed/Root traits/Exp2/Coord data for cicling figure.txt/", header = TRUE)
str(Coordata)
Coordata$Microbiome = factor(Coordata$Microbiome, c('SynCom', 'NB'))
Coordata$Cell_type = factor(Coordata$Cell_type, c('Epidermis','Exodermis','CCL4','CCL3','CCL2','CCL1','Endodermis','Stele')) 
Coordata$Species = factor(Coordata$Species, c('LY5280','SP9509','SI7820','LP0049','LJ9250','LA7339','LP8539')) 

Coordata2=Coordata[Coordata$Experiment=='1',]
CoordataV2 = Coordata[Coordata$Nutrition!='1.34mM',]
CoordataV3 = CoordataV2[CoordataV2$Species!='LY5280',]

Mean.Traits <- CoordataV3%>%group_by(Nutrition,Microbiome,Species,Cell_type)%>%
  summarise_all(mean)

SD.Traits <- CoordataV3[,-c(3,5)]%>%group_by(Nutrition,Microbiome,Species,Cell_type)%>%
  summarise_all(se)

colnames(SD.Traits) = c("Nutrition" ,"Microbiome","Species", "Cell_type","Subsample","Experiment","Absoluto_se","Percetange_se")
Mean.Traits = data.frame(cbind(Mean.Traits, SD.Traits[,7:8])) 


Mean.TraitsB = Mean.Traits[Mean.Traits$Microbiome!='NB',]
Mean.TraitsBHigh = Mean.TraitsB[Mean.TraitsB$Nutrition=='335mM',]
Mean.TraitsBMedium = Mean.TraitsB[Mean.TraitsB$Nutrition=='6.7mM',]
Mean.TraitsBLow = Mean.TraitsB[Mean.TraitsB$Nutrition=='0.67mM',]
Mean.TraitsBLow[Mean.TraitsBLow$Species=='LP8539',]

Mean.TraitsNB = Mean.Traits[Mean.Traits$Microbiome=='NB',]
Mean.TraitsNBHigh = Mean.TraitsNB[Mean.TraitsNB$Nutrition=='335mM',]
Mean.TraitsNBMedium = Mean.TraitsNB[Mean.TraitsNB$Nutrition=='6.7mM',]
Mean.TraitsNBLow = Mean.TraitsNB[Mean.TraitsNB$Nutrition=='0.67mM',]
Mean.TraitsNBLow[Mean.TraitsNBLow$Species=='LP8539',]

xpBHigh = ggplot(data=Mean.TraitsBHigh, aes(x=Species, y=Percentage, fill=Cell_type)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme_bw()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") + 
  ylim(0,101)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))

pBCircHigh = pBHigh + coord_polar()+ theme(legend.position="none")+xlab('')

H_B = ggplot(data=Mean.TraitsBHigh, aes(x=Species, y=Percentage, fill=Cell_type)) + #pBCirc =
  theme(axis.text.x = element_text(angle = 45, hjust=1,),
        legend.position="left",
        strip.background =  element_rect(colour = "black"))+
  theme_linedraw()+
  geom_bar(stat="identity", color = "black") +
  #coord_polar()+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))+
  geom_errorbar(aes(ymin=Percentage, ymax=Percentage+Percetange_se), 
                position=position_dodge(0.5), width=1.5, alpha = 1, stat = 'identity') +
  labs(title = "Microbiome",y='Cell diameter (%)',x='Species')

pBCircMedium = pBMedium + coord_polar()+ theme(legend.position="none")+xlab('')

pBMedium = ggplot(data=Mean.TraitsBMedium, aes(x=Species, y=Percentage, fill=Cell_type)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme_bw()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,101)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))

M_B = ggplot(data=Mean.TraitsBMedium, aes(x=Species, y=Percentage, fill=Cell_type)) + #pBCirc =
  theme(axis.text.x = element_text(angle = 45, hjust=1,),
        legend.position="left",
        strip.background =  element_rect(colour = "black"))+
  theme_linedraw()+
  geom_bar(stat="identity", color = "black") +
  #coord_polar()+
  scale_fill_manual(values=c("#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))+
  geom_errorbar(aes(ymin=Percentage, ymax=Percentage+Percetange_se), 
                position=position_dodge(0.5), width=1.5, alpha = 1, stat = 'identity') +
  labs(title = "Microbiome",y='Cell diameter (%)',x='Species')

pBLow = ggplot(data=Mean.TraitsBLow, aes(x=Species, y=Percentage, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme_bw()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,101)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pBCircLow = pBLow + coord_polar()+ theme(legend.position="none")+xlab('')

L_B = ggplot(data=Mean.TraitsBLow, aes(x=Species, y=Percentage, fill=Cell_type)) + #pBCirc =
  theme(axis.text.x = element_text(angle = 45, hjust=1,),
        legend.position="left",
        strip.background =  element_rect(colour = "black"))+
  theme_linedraw()+
  geom_bar(stat="identity", color = "black") +
  #coord_polar()+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))+
  geom_errorbar(aes(ymin=Percentage, ymax=Percentage+Percetange_se), 
                position=position_dodge(0.5), width=1.5, alpha = 1, stat = 'identity') +
  labs(title = "Microbiome",y='Cell diameter (%)',x='Species')

pNBHigh = ggplot(data=Mean.TraitsNBHigh, aes(x=Species, y=Percentage, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme_bw()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,101)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pNBCircHigh = pNBHigh + coord_polar()+ theme(legend.position="none")+xlab('')

H_nB = ggplot(data=Mean.TraitsNBHigh, aes(x=Species, y=Percentage, fill=Cell_type)) + #pBCirc =
  theme(axis.text.x = element_text(angle = 45, hjust=1,),
        legend.position="left",
        strip.background =  element_rect(colour = "black"))+
  theme_linedraw()+
  geom_bar(stat="identity", color = "black") +
  #coord_polar()+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))+
  geom_errorbar(aes(ymin=Percentage, ymax=Percentage+Percetange_se), 
                position=position_dodge(0.5), width=1.5, alpha = 1, stat = 'identity') +
  labs(title = "Microbiome",y='Cell diameter (%)',x='Species')

pNBMedium = ggplot(data=Mean.TraitsNBMedium, aes(x=Species, y=Percentage, fill=Cell_type)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme_bw()+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") +ylim(0,101)+
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pNBCircMedium = pNBMedium + coord_polar()+ theme(legend.position="none")+xlab('')

M_nB = ggplot(data=Mean.TraitsNBMedium, aes(x=Species, y=Percentage, fill=Cell_type)) + #pBCirc =
  theme(axis.text.x = element_text(angle = 45, hjust=1,),
        legend.position="left",
        strip.background =  element_rect(colour = "black"))+
  theme_linedraw()+
  geom_bar(stat="identity", color = "black") +
  #coord_polar()+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))+
  geom_errorbar(aes(ymin=Percentage, ymax=Percentage+Percetange_se), 
                position=position_dodge(0.5), width=1.5, alpha = 1, stat = 'identity') +
  labs(title = "Microbiome",y='Cell diameter (%)',x='Species')

pNBLow = ggplot(data=Mean.TraitsNBLow, aes(x=Species, y=Percentage, fill=Cell_type)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme_bw(aes(color='black'))+ 
  theme(legend.position="left") +
  geom_bar(stat="identity", color = "black") + ylim(0,101) + 
  scale_fill_manual(values=c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))
pNBCircLow = pNBLow + coord_polar()+ theme(legend.position="none")+xlab('')

L_nB = ggplot(data=Mean.TraitsNBLow, aes(x=Species, y=Percentage, fill=Cell_type)) + #pBCirc =
  theme(axis.text.x = element_text(angle = 45, hjust=1,),
        legend.position="left",
        strip.background =  element_rect(colour = "black"))+
  theme_linedraw()+
  geom_bar(stat="identity", color = "black") +
  #coord_polar()+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))+
  geom_errorbar(aes(ymin=Percentage, ymax=Percentage+Percetange_se), 
                position=position_dodge(0.5), width=1.5, alpha = 1, stat = 'identity') +
  labs(title = "Microbiome",y='Cell diameter (%)',x='Species')

ggpubr::ggarrange(pNBCircHigh,pBCircHigh,pNBCircMedium,pBCircMedium,pNBCircLow,pBCircLow, ncol = 2, nrow=3)

ggpubr::ggarrange(H_nB,H_B,M_nB,M_B,L_nB,L_B, ncol = 2, nrow=3,common.legend = TRUE, legend = 'left')

#circle error

