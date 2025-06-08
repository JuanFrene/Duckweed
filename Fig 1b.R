########## Bar plot de las diferencias de  cell layers
#Tables
diff_cell_layer <- read.table("clipboard", header = TRUE)
diff_cell_layer$Species <-factor(diff_cell_layer$Species,c('PS','LY5280','SI9227','SP7498','SP9509',
                                                           'SP7820','SP5543','SP9192','LM7200','LT9243',
                                                           'LP7760','LV7650','LJ9250','LP0049','LM8389','LT9109',
                                                           'LA7339','LM7123','LP8539')) 

ggplot(diff_cell_layer, aes(Species, diff_cell_layer,  fill = Complexity)) +
  geom_boxplot() + 
  theme_few() +
  geom_hline(yintercept = 0) +
  labs(y="Difference of cell layers", 
       x="Species") 

#Radar figure
Root_Trait_Mean <- TablaCompleta%>%group_by(Species, Microbiome)%>%
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
