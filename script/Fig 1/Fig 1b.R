setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/Exp1/")

Root_traits = read.table("Root traits.txt", header = TRUE, row.names = 1)

Root_Trait_Mean <- Root_traits%>%group_by(Species, Microbiome)%>%
  summarise_all(mean)
Root_Trait_Mean2=data.frame(Root_Trait_Mean)

ggplot(data=Root_Trait_Mean2,  aes(x=Species, y=Still_diameter , group=Microbiome, colour=Microbiome)) + 
  geom_point() + 
  geom_line() + 
  xlab("Species") + 
  ylab("um") + 
  ylim(0,45) + ggtitle("Steele Diameter")  + 
  geom_hline(aes(yintercept=0), lwd=1, lty=2) +
  coord_polar()+
  theme_bw()

ggplot(data=Root_Trait_Mean2[Root_Trait_Mean2$Species!="PS",],  aes(x=Species, y=Endodermis_diameter , group=Microbiome, colour=Microbiome)) + 
  geom_point() + 
  geom_line() + 
  xlab("Species") + 
  ylab("um") + 
  ylim(0,10) + ggtitle("Stele Diameter")  + 
  geom_hline(aes(yintercept=0), lwd=1, lty=2) +
  coord_polar()+
  theme_bw()

ggplot(data=Root_Trait_Mean2,  aes(x=Species, y=Arenquima_Perc, group=Microbiome, colour=Microbiome)) + #[Root_Trait_Mean2$Species!="PS",]
  geom_point() + 
  geom_line() + 
  xlab("Species") + 
  ylab("%") + 
  ylim(0,30) + ggtitle("Aerenchyma Perc. Area")  + 
  geom_hline(aes(yintercept=0), lwd=1, lty=2) +
  coord_polar()+
  theme_bw()

