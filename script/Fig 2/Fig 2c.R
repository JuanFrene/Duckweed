packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'reshape2')
sapply(packages, require, character.only = TRUE)              

Alphadiversity = read.table("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 2/Fig 2c/Alphadiversity en txt.txt", header = TRUE, row.names = 1)
Alphadiversity2 = Alphadiversity[Alphadiversity$Shannon > 0.01,]

shannon.model.Root<-lmer(Shannon ~ Specie  + (1|Sample), data = Alphadiversity2[Alphadiversity$Root=='Low',])
aov = anova(shannon.model.Root)

#High
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Specie 0.35904 0.039893     9 23.555  2.5697 0.03196 *

#Medium
#Type III Analysis of Variance Table with Satterthwaite's method
#       Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)   
#Specie 0.57581 0.063979     9 24.568  4.5023 0.00142

#Low
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Specie 0.91438  0.1016     9 29.131  5.3306 0.0002534 ***
---

#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF DenDF  F value    Pr(>F)    
#Root     0.99539 0.49769     2 106.7  19.313    6.941e-08


paleta_alive <- c("#C200FF",'#FF0000','#8B7500')
Alphadiversity2$Root = factor(Alphadiversity2$Root, c('High', 'Medium', 'Low'))

ggplot(Alphadiversity2, aes(Species, Shannon), ) + 
  geom_boxplot(aes(fill=Complexity)) +
  facet_grid(.~Root)+
  theme_few()+
  scale_fill_gradientn(colours=c('blue','red'))+ 
  labs(title ="Shannon index", x='Synthetic Root') + 
  theme(legend.position="right",axis.text.x = element_text(angle=65, vjust=0.6))
  ylim(0,NA)
  
ggplot(data = Alphadiversity2, aes(Complexity, Shannon)) +
    #geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.2, alpha = 0.4) +
    #geom_errorbar(mapping = aes(ymin =mean.x - ds.x ,ymax = mean.x + ds.x),size = 0.2, alpha = 0.4)  + 
    #geom_hline(yintercept = 0,size = 0.2,color = "black",linetype = "longdash")+
  facet_grid(.~Root)+
  geom_point( size = 2, aes(colour = Specie),stroke = 1) + 
    geom_smooth(method="lm",size = 0.5, se=T,,color='black') + 
    labs(y='PCoA 2', x='CCL') +
    theme_few()+
    scale_color_manual(values = paleta_alive)
  
summary(lm(Shannon~Complexity, data = Alphadiversity2))

  