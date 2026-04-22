library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(harrietr)
library(egg)
library(paletteer)
library(scales)
library(car)
library(Rmisc)
library(reshape2)
library(ggthemes)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp1/")

Status <- read.table("Status2.txt", header = TRUE, row.names = 1)

Status$Treatment <- factor(Status$Treatment, c("With-Microbiome","No-Microbiome"))
Status$Species <- factor(Status$Species,rev(c('PS','SI7820','SI9227','SP9509','SP5543',
                                          'SP7498','SP9192','LP7760','LP0049',
                                          'LV7650','LJ9250','LM7123','LM8389','LT9109',
                                          'LY5280','LM7200','LA7339','LP8539','LT9243')))

head(Status)
ggplot(data=Status, aes(y=Species, x=Dead, fill=Treatment )) +
  geom_boxplot(outlier.shape = NA)+
  #geom_bar(stat="identity", position=position_dodge(), color = 'black', size = 0.3,width = 0.9) +# 
  scale_fill_manual(values = c('orange','#0000ff')) +
  theme_few() + guides() + #color = FALSE, fill=FALSE
  labs(x="Mortality (%)", y = 'Species') +
  theme(legend.title=element_blank(), 
        legend.margin=margin(c(0,0,0,0)))

  