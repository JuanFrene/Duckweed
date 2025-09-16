setwd("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig S1/Fig S1b/")

Data <- read.table("number of roots.txt",header = TRUE , row.names = 1)

Data$Species <-factor(Data$Species, c('LT9243','LP8539','LA7339','LM7200','LY5280','LV7650',
                                      'LM8389','LM7123','LJ9250','LT9109','LP0049','LP7760',
                                      'SP9192','SP7498','SP5543','SP9509','SI9227','SP7820','PS'))


ggplot(Data, aes(Species, root_plant,  fill = Treatment)) + 
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_few()+ 
  labs(y='Number of roots',x='Species')+
  scale_fill_manual(values=c('orange', '#0000ff'),l=c('Sterile Water','Natural Water'))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  theme(legend.position="top")+
  ylim(0,10)
