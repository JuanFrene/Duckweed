packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')
sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/") #Users/juanp/Google

seqtab_Root <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 2/Fig 2a/seqtab_Root.rds")

###CAP
CAP2 <- ordinate(seqtab_Root, "CAP", "bray",~data$Nutrient)
summary(CAP2)
xxlab_text = 'CAP1 (46.63%)'
ylab_text = 'CAP2 (23.86%)'

x<-data.frame(CAP2$CCA$wa[, 1:2])
x <- merge(x,seqtab_Root@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
nrow(x)

safe_colorblind_palette <- c("#44AA99", "#999933", "#882255", "#332288")

safe_colorblind_palette2 <- c("#1912f1", "#117733","#DDCC77","#882615" )

x$Complexity <- "High"
x$Complexity[which(x$Species == 'LA7339')] <- "Low"
x$Complexity[which(x$Species == 'LP8539')] <- "Low"
x$Complexity[which(x$Species == 'LJ9250')] <- "Medium"
x$Complexity[which(x$Species == 'LP0049')] <- "Medium"
x$Complexity <- x$Complexity %>% factor

ggplot(x, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Nutrient), size=4,pch=21, color="black") +
  theme_few() +
  geom_tile()+
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  #labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette)+
  geom_text(aes(label = "p = 0.776",
    x = -0.6, y = 0.63, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.032",
                x = -0.583, y = 0.57, size = 6, fill = "black"))

ggplot(x, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Complexity ), size=4,pch=21, color="black") +
  theme_few() +
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  #labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette2)+
  geom_text(aes(label = "p = 0.049",
                x = -0.6, y = 0.63, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.038",
                x = -0.587, y = 0.57, size = 6, fill = "black"))

