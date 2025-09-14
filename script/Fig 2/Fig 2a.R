packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')
sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/") #Users/juanp/Google

seqtab_Root <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 2/Fig 2a/seqtab_Root.rds")
seqtab_Frond <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 2/Fig 2a/seqtab_Frond.rds")

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



###CAP Frond
dataF = seqtab_Frond@sam_data
CAPF <- ordinate(seqtab_Frond, "CAP", "bray",~dataF$Nutrient)
summary(CAPF)
xlab_text = 'CAP1 (33,06%)'
ylab_text = 'CAP2 (11.49%)'

f<-data.frame(CAPF$CCA$wa[, 1:2])
f <- merge(f,seqtab_Frond@sam_data, by=0)
rownames(f)<-x$Row.names
f$Row.names <- NULL
nrow(f)

f$Complexity <- "High"
f$Complexity[which(f$Species == 'LA7339')] <- "Low"
f$Complexity[which(f$Species == 'LP8539')] <- "Low"
f$Complexity[which(f$Species == 'LJ9250')] <- "Medium"
f$Complexity[which(f$Species == 'LP0049')] <- "Medium"
f$Complexity <- f$Complexity %>% factor

safe_colorblind_palette <- c("#44AA99", "#999933", "#882615", "#332288", "#6699CC", "#888888")

ggplot(f, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Nutrient), size=4,pch=21, color="black") +
  theme_few() +
  geom_tile()+
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette)+
  geom_text(aes(label = "p = 0.008",
                x = -0.443, y = 0.66, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.044",
                x = -0.43, y = 0.6, size = 6, fill = "black"))

ggplot(f, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Complexity ), size=4,pch=21, color="black") +
  theme_few() +
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette2)+
  geom_text(aes(label = "p = 0.297",
                x = -0.443, y = 0.66, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.013",
                x = -0.433, y = 0.6, size = 6, fill = "black"))



