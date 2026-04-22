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
x <- merge(x, seqtab_L@sam_data, by=0)
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

library(paletteer)
library(pals)


  ###PERMANOVA
# Calculate bray curtis distance matrix
ps.3_bray <- phyloseq::distance(seqtab_C, method = "bray")
# make a data frame from the sample_data
sample_ps.3c <- data.frame(sample_data(seqtab_C))
# Adonis test
PERMANOVA_C = adonis2(ps.3_bray ~ Root*Species , data = sample_ps.3c, by = "terms" )

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
PermanovaL = adonis2(ps.Lysine_bray ~ Root*Species , data = sample_ps.Lysine)

#Lysine: 
#            Df SumOfSqs      R2      F Pr(>F)   
#Root          2   0.4731 0.08792 1.0734  0.359   
#Species       1   0.8163 0.15170 3.7041  0.003 **
#Root:Species  2   0.3453 0.06416 0.7834  0.662   
#Residual     17   3.7465 0.69622                 
#Total        22   5.3813 1.00000   

#PERMANOVA
set.seed(999)

PERMANOVA2 = data.frame(PERMANOVA_C)
variable  = c('variable','variable','variable','variable','variable')
PERMANOVA2$Significance <- "No Significant"
pval_thres <- 0.05
PERMANOVA2$Significance[which(PERMANOVA2$Pr..F. < pval_thres)] <- "Significant"
PERMANOVA2$Significance <- PERMANOVA2$Significance %>% factor

LC = data.frame(cbind(variable, row.names(PERMANOVA2), PERMANOVA2[,c(3,6)]))
colnames(LC) = c('variable','Effect', 'R2','Significance')
Lysi = data.frame(cbind(variable  = c('variable','variable','variable','variable'),
                Effect = c('Root','Species','Root:Species','Residual'),
                R2 = c(0.08792,0.15170,0.06416,0.69622),
                Significance = c("No Significant","Significant","No Significant","No Significant")))

Lysi$R2 = as.numeric(Lysi$R2)

Cont = data.frame(cbind(variable  = c('variable','variable','variable','variable'),
                        Effect = c('Root','Species','Root:Species','Residual'),
                        R2 = c(0.25706,0.02858,0.03870,0.67565),
                        Significance = c("Significant","No Significant","No Significant","No Significant")))

Cont$R2 = as.numeric(Cont$R2)

#Plot Variance
dim(LC)
colnames(LC)
dim(LC)
data_melt <- Cont

#Reorder x axis: the order of the factor will be the same as in the data.csv file

data_melt$Effect <- as.character(data_melt$Effect)#Turn your 'treatment' column into a character vector
data_melt$Effect <- factor(data_melt$Effect, levels=unique(rev(data_melt$Effect)))#Then turn it back into a factor with the levels in the correct order

mypal = c('white','#c0047b','#a9d11c','#1f4f58','#528501')

#data_melt$Effect = factor(data_melt$Effect, c('Crop','Compost','Crop:Compost','Residual'))

ggplot(data=data_melt, aes(x=variable, y=R2,fill=Effect)) +
  geom_bar(stat="identity",aes(color=Significance), size = 0.3,width = 1,linewidth=1) +# 
  scale_fill_manual(values = mypal) +
  theme_few() + #guides() + #color = FALSE, fill=FALSE
  labs(y="Variance explained (%)") +
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significance vs Water")+
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0)),
        axis.text.x   = element_blank())

