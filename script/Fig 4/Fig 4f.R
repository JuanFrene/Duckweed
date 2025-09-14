packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')

sapply(packages, require, character.only = TRUE)              

seqtab_C <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 4/Fig 4f/Control.rds")
seqtab_L <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 4/Fig 4f/Lysine.rds")

#####PCoA for ASV-level data with Bray-Curtis
PCoA <- ordinate(seqtab_C, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x, seqtab_C@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Species*Root,x,mean)
PCoAmeands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Species*Root,x,se)
PCoAmeanLJ =merge(PCoAmean[c(1,4),],PCoAmeands[c(1,4),-(1)],by="Root")
PCoAmeanSP =merge(PCoAmean[-c(1,4),],PCoAmeands[-c(1,4),-(1)],by="Root")
PCoAmean2 = rbind(PCoAmeanLJ,PCoAmeanSP)


paleta_alive <- c("#C200FF",'#8B7500','#FF0000','#00008B',"#FFB919",'#FF7F50',"#00CC7A")

PCoAmean2$Root = factor(PCoAmean2$Root, c('High', 'Medium','Low'))
ggplot(data = PCoAmean2, aes(mean.x, mean.y, label=Species)) +
    geom_vline(xintercept = 0, size = 2, color = "#D9D9D9",linetype = "longdash")+
    geom_hline(yintercept = 0, size = 2, color = "#D9D9D9",linetype = "longdash")+
    geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),linewidth = 0.1, alpha = 0.15) +
    geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),linewidth = 0.1, alpha = 0.15)  + 
    geom_point( size = 6, aes(colour = Root, shape = Species),stroke = 1) + 
    #scale_fill_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
    labs(y=ylab_text, x=xlab_text)+ theme_few()+
    scale_colour_manual(values = paleta_alive)

ggplot(x, aes(Axis.1, Axis.2)) + #[PCoAmean$Specie!='PS',]
  geom_point(aes(color=Species,shape=Root, size=6)) +
  #geom_text(hjust=0.5, vjust=-0.5) + 
  theme_few() +
  theme(legend.position="top")+ 
  scale_fill_paletteer_d("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  #scale_colour_gradientn(colours=c("Blue","Red"))+
  labs(title = "Bacterial community associated with Duckweed Root",
       y=ylab_text, 
       x=xlab_text)    

library(paletteer)
library(pals)


  ###PERMANOVA
# Calculate bray curtis distance matrix
ps.3_bray <- phyloseq::distance(seqtab_C, method = "bray")
# make a data frame from the sample_data
sample_ps.3c <- data.frame(sample_data(seqtab_C))
# Adonis test
adonis2(ps.3_bray ~ Root*Species , data = sample_ps.3c)

#Control            Df SumOfSqs      R2      F Pr(>F)  
#Root          2   1.1450 0.25706 2.6632  0.020 *
#Species       1   0.1273 0.02858 0.5923  0.677  
#Root:Species  1   0.1724 0.03870 0.8019  0.519  
#Residual     14   3.0095 0.67565                
#Total        18   4.4541 1.00000 


ps.Lysine_bray <- phyloseq::distance(seqtab_L, method = "bray")
# make a data frame from the sample_data
sample_ps.Lysine <- data.frame(sample_data(seqtab_L))
# Adonis test
adonis2(ps.Lysine_bray ~ Root*Species , data = sample_ps.Lysine)

#Lysine: 
#            Df SumOfSqs      R2      F Pr(>F)   
#Root          2   0.4731 0.08792 1.0734  0.359   
#Species       1   0.8163 0.15170 3.7041  0.003 **
#Root:Species  2   0.3453 0.06416 0.7834  0.662   
#Residual     17   3.7465 0.69622                 
#Total        22   5.3813 1.00000   

