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

AlfaRootMean = aggregate(Shannon~Specie.1,AlfaCorrb3[AlfaCorrb3$Compartment=='Root',],mean)
AlfaFrontMean = aggregate(Shannon~Specie.1,AlfaCorrb3[AlfaCorrb3$Compartment=='Front',],mean)

mortality = c(68.85245902,9.677419355,14.94252874,34.375,21.05263158,8.421052632,14.45783133,63.51351351,63.38028169,39.63963964,45,5.376344086,0,1.315789474,1.282051282,12,7.291666667,4.651162791,15.55555556)
RootComp = c(5,5.33,5,6.4,5,5.2,6,5,5.2,6.2,6,7.6,8.4,7.6,6.6,7.25,6.8,6.6,7)
RootComp2 = c(5,5,5.75,6.75,5,5,5.2,5.25,5,5.4,5.6,5.8,7.6,7.333333333,6,6.8,5.8,6,6.6)
Epidermis = c(16.77,14.41,15.3,10.44,18.79,8.07,10.05,14.31,17.55,12.54,7,11.49,32.84,8.06,10.89,10.45,10.44,8.37,9.17)
EpidermisPorc = c(51.86,23.27,17.73,22.56,19.92,23.66,48.81,21.44,26.83,27.71,26.04,18.71,32.88,24.19,41.42,40.69,38.31,39.09,32.25)
Nroot = c(2,1.478333333,1.285714286,1,1,1.946666667,2.266666667,1,1,1,1,1,5.8,2.893333333,3.12,2.406666667,2.82,2.341666667,3.473333333)


AlfaRootMean = cbind(AlfaRootMean, mortality,RootComp,Epidermis,Nroot)
species = data.frame(AlfaRootMean[,1])
colnames(species) = 'Specie'
AlfaRootMean2 = cbind(species, RootComp, mortality,Nroot)
AlfaFrontMean = cbind(AlfaFrontMean, mortality,RootComp,Epidermis,Nroot)

CorralfaDeadR <- ggplot(AlfaRootMean, aes(RootComp, Shannon)) + 
  labs(title="Root alphadiversity", y="Shannon index ", x= 'Root Complexity') +
  ylim(0,5)

CorralfaDeadR + 
  geom_jitter(aes(color=mortality, size=6)) + 
  theme_few() +
  geom_smooth(method="lm", se=T)  + 
  scale_colour_gradientn(colours=c("Dark Green","yellow"))+
  geom_text(aes(x = 5.5, y = 5,
                label = "R2 = 0.2599, P = 0.02577"),
            stat = "unique")

CorralfaFroxntCompF <- ggplot(AlfaFrontMean, aes(RootComp, Shannon)) + 
  labs(title="Frond alphadiversity", y="Shannon index ", x= 'Root Complexity') +
  ylim(0,5)

CorralfaFroxntCompF + 
  geom_jitter(aes(color=mortality, size=6 )) + 
  theme_few() +
  geom_smooth(method="lm", se=T) + 
  scale_colour_gradientn(colours=c("Dark Green","yellow"))+
  geom_text(aes(x = 5.5, y = 5,
                label = "R2 = 0.0591, P = 0.3465"),
            stat = "unique")
