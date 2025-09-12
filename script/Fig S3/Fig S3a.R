packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest",'dada2', "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", "reshape2", 'lme4')
sapply(packages, require, character.only = TRUE)              
setwd("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/")

Alphadiver <- read.table("Alphadiv data.txt", header = TRUE)
Metadata <- read.table("Metadata.txt", header = TRUE)

AlfaCorr1 <- merge('Metadata','Alphadiver', by='nowmanes')

AlfaCorrb2$Specie.1   <- factor(AlfaCorrb2$Specie.1    ,c('PS','LY5280','SI9227','SP7498','SP9509',
                                                          'SI7820','SP5543','SP9192','LM7200','LT9243',
                                                          'LP7760','LV7650','LJ9250','LP0049','LM8389','LT9109',
                                                          'LA7339','LM7123','LP8539', 'Control', 'Inoculum'))

AlfaCorr1$Compartment   <- factor(AlfaCorr1$Compartment, c('Inoculum','Front', 'Root','Water'))

AlfaCorr2 <- AlfaCorr1[AlfaCorr1$Compartment != 'Inoculum',]
AlfaCorr3 <- AlfaCorr2[AlfaCorr2$Specie.1 != 'Control',]

summary(shannon.Compartment.root<-aov(Shannon ~Specie.1, data = AlfaCorrb[AlfaCorrb$Compartment=='Front',]))
lsdshannon.Compartment.root <- HSD.test(shannon.Compartment.root,'Specie.1', alpha = 0.05, group=TRUE) #, p.adj="bonferroni"
lsdshannon.Compartment.root


AlfaCorrb= data.frame(AlfaCorr3[AlfaCorr3$Shannon>0.3,])  
AlfaCorrb$Shannon = as.numeric(AlfaCorrb$Shannon)


AlfaCorrb$Specie.1   <- factor(AlfaCorrb$Specie.1    ,c('LT9243','LP8539','LA7339','LM7200','LY5280','LV7650',
                                                        'LM8389','LM7123','LJ9250','LT9109','LP0049','LP7760',
                                                        'SP9192','SP7498','SP5543','SP9509','SI9227','SI7820','PS'))

ggplot(AlfaCorrb[AlfaCorrb$Compartment!='Water',], aes(Shannon, Specie.1)) + 
  geom_boxplot(aes(fill=Compartment)) +
  facet_grid(.~Compartment)+
  theme_few()+ scale_fill_manual(values=c('#FFA000', '#458CFF'))+
  theme(legend.position="top")
