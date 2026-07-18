packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')
sapply(packages, require, character.only = TRUE)              


seqtab_Frond <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig S16/Fig S16a/seqtab_Frond.rds")

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



