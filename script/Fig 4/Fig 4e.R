packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')

sapply(packages, require, character.only = TRUE)              

seqtab <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig 4/Fig 4e/ps.perc.Lysine.rds")

#####PCoA for ASV-level data with Bray-Curtis
# Calculates Bray-Curtis distances between samples. Because taxa is in
# columns, it is used to compare different samples. We transpose the
# assay to get taxa to columns
PCoA <- ordinate(seqtab, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x, seqtab@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Species*Drugs,x,mean)
PCoAmeands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Species*Drugs,x,se)
PCoAmean2 =merge(PCoAmean,PCoAmeands[,-(2)],by="Species")

gg2 <- merge(PCoAmean2,x[,c(9,12)], by="Species")
gg2$Complexity = as.numeric(gg2$Complexity)

ggplot(data = gg2, aes(mean.x, mean.y, label=Species)) +
    geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
    geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
    geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.1, alpha = 0.15) +
    geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),size = 0.1, alpha = 0.15)  + 
    geom_point( size = 6, aes(col = Complexity),stroke = 1) + 
    scale_colour_gradientn(colours=c("Blue","Red"))+
    labs(y=ylab_text, x=xlab_text)+ theme_few()
    
###PERMANOVA
# Calculate bray curtis distance matrix
ps.3_bray <- phyloseq::distance(seqtab, method = "bray")
# make a data frame from the sample_data
sample_ps.3c <- data.frame(sample_data(seqtab))
# Adonis test
adonis2(ps.3_bray ~ Species, data = sample_ps.3c)


#Lysine: 
#adonis2(formula = ps.3_bray ~ Complexity * Species, data = sample_ps.3c)
#Complexity  1   0.5616 0.14202 2.2936  0.060 .
#Species     2   0.6995 0.17688 1.4283  0.144  
#Residual   11   2.6936 0.68110                

