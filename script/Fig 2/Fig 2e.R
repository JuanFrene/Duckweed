packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'reshape2')
sapply(packages, require, character.only = TRUE)              
library(dada2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Microtubulos/16S/")

###ASV table
asv.table2 <- readRDS('G:/My Drive/labs/Nottingham/Duckweed/Exp Microtubulos/16S/seqtab_final.rds')
asv.table<- otu_table(asv.table2, taxa_are_rows=FALSE)

##ASing to the syncom
taxa <- assignTaxonomy(asv.table, "G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/Syncom_Sequence.pcr.unique.fasta", multithread=TRUE)#, minBoot =97

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

#Now we can make the phyloseq object
ps <- phyloseq(asv.table, tax_table(taxa), sample_data(meta))
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps.1 = subset_taxa(ps, Genus  != "NA")

ps.2 = subset_samples(ps.1, ID !=  "Blanck")
ps.2.pruned <- prune_taxa(taxa_sums(ps.2)>=1, ps.2)
ps.2.pruned2 <- prune_samples(sample_sums(ps.2.pruned)>0, ps.2.pruned) 

ps.3 = subset_samples(ps.2, Root !=  "OR")
ps.3.pruned <- prune_taxa(taxa_sums(ps.3)>=1, ps.3)
ps.3.pruned2 <- prune_samples(sample_sums(ps.3.pruned)>0, ps.3.pruned) 

ps.4 = subset_samples(ps.3, ID !=  "Syncom")
ps.4.pruned <- prune_taxa(taxa_sums(ps.4)>=1, ps.4)
ps.4.pruned2 <- prune_samples(sample_sums(ps.4.pruned)>0, ps.4.pruned) 

ps.SR2 = subset_samples(ps.4.pruned2, Compartment ==  "SR")
ps.SR1 <- prune_taxa(taxa_sums(ps.SR2)>=1, ps.SR2)
ps.SR <- prune_samples(sample_sums(ps.SR1)>0, ps.SR1) 
ps.SR.H <- transform(ps.SR, "hellinger")
ps.SR.perc <- transform_sample_counts(ps.SR, function(x) x / sum(x)) 

ps.M2 = subset_samples(ps.4.pruned2, Compartment !=  "SR")
ps.M1 <- prune_taxa(taxa_sums(ps.M2)>=1, ps.M2)
ps.M <- prune_samples(sample_sums(ps.M1)>0, ps.M1) 


####PCoA
ps.SR.Low = subset_samples(ps.SR, Root ==  "Low")
ps.SR.Medium = subset_samples(ps.SR, Root ==  "Medium")
ps.SR.High = subset_samples(ps.SR, Root ==  "High")

PCoA <- ordinate(ps.SR.High, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100

xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x, ps.SR@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
x$ID


PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~ID,x,mean) #
PCoAmeands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~ID,x,se)
PCoAmean2 =merge(PCoAmean,PCoAmeands,by="ID")
PCoAmean2[is.na(PCoAmean2)] <- 0.03
PCoAmean3 =merge(PCoAmean2,meta[,-(5)],by="ID")
PCoAmean4 <- PCoAmean3%>%group_by(ID, Species, Compartment, Root)%>%
  summarise_all(mean)


paleta_alive <- c("#C200FF",'#FF0000','#8B7500','#8B8B7A',"#FFB919",'#FF7F50',"#00CC7A",'#8B1A1A','#104E8B','#8B1C62','#7FFF00')


###Control alone
x$Species = factor(x$Species, c('LY9205','SI9227','SP7498','SP9509','SI7820','LT9243','LJ9250','LP0049','LA7339','LP8539','Control'))
x$Root = factor(x$Root, c('High','Medium','Low'))

ggplot(data = x[x$Species=='Control',], aes(Axis.1, Axis.2)) + #[x$Species=='Control',]
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(size = 2, aes(color = Root),stroke = 1) + #, shape=Compartment
  #facet_grid(.~Species, space = "fixed",scales = "fixed") +
  xlab(label = xlab_text) + ylab(label = ylab_text) +
  scale_fill_manual(values = paleta_alive) +
  scale_shape_manual(values=c(0, 1, 2,3,4,5,6,7,8,9,10))+
  theme_few()+
  scale_color_manual(values = paleta_alive)+
  theme(plot.background = element_rect(size = 2))+ 
  theme(panel.border = element_rect( size=2))


PCoAmean4$Species = factor(PCoAmean4$Species, c('LY9205','SI9227','SP7498','SP9509','SI7820','LT9243','LJ9250','LA7339','LP0049','LP8539','Control'))
PCoAmean4$Root = factor(PCoAmean4$Root, c('High','Medium','Low'))
PCoAmean4 = data.frame(PCoAmean4)
PCoAmean4$Complexity[is.na(PCoAmean4$Complexity)] =  6

PCoAmean4$Complexity = as.numeric(PCoAmean4$Complexity)

ggplot(data = PCoAmean4[PCoAmean4$Species!="Control",], aes(mean.x, mean.y, name=Species)) +theme_few()+
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash") +
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash") +
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),size = 0.02,alpha = 0.2) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),size = 0.02,alpha = 0.2)  + 
  facet_grid(.~Root, space = "fixed",scales = "fixed") +
  geom_point(aes(colour = Complexity),stroke = 1,shape = 19,size=3) + #, size=3,
  xlab(label = xlab_text) + ylab(label = ylab_text) +
  scale_colour_gradientn(colours=c("Blue","Red"))+
  scale_shape_manual(values=c(0, 1, 2,3,4,5,6,7,8,9,10))+
  theme(panel.border = element_rect( size=2))


###CAP
ps.SR.sc = subset_samples(ps.SR, Species !=  "Control")
data = ps.SR.sc@sam_data

CAP2 <- ordinate(ps.SR.sc, "CAP", "bray",~data$Root)
summary(CAP2)
xxlab_text = 'CAP1 (46.63%)'
ylab_text = 'CAP2 (23.86%)'

x<-data.frame(CAP2$CCA$wa[, 1:2])
x <- merge(x,data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
nrow(x)

ggplot(x, aes(CAP1, CAP2)) + #
  geom_point(aes(colour = Root, shape = Root), size=4) +
  theme_few() +
  geom_tile()+
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = paleta_alive)

###PERMANOVA
ps.SR.Low = subset_samples(ps.SR.sc, Root ==  "Low")
ps.SR.Medium = subset_samples(ps.SR.sc, Root ==  "Medium")
ps.SR.High = subset_samples(ps.SR.sc, Root ==  "High")

# Calculate bray curtis distance matrix
ps.3_bray2 <- phyloseq::distance(ps.SR.Low, method = "bray")
# make a data frame from the sample_data
sample <- data.frame(sample_data(ps.SR.Low))
# Adonis test
PERMANOVA = adonis2(ps.3_bray2 ~ Complexity*Species , data = sample)

#PERMANOVA_H_C Complexity:   R2 = 0.00606 F=0.2011  P=0.991
#PERMANOVA_H_S Species:    R2 = 0.30234 F=1.4084  P=0.092 .
#PERMANOVA_M_C Complexity: R2 = 0.03353 F=1.1797  P=0.305
#PERMANOVA_M_S Species: R2 = 0.21444 F=0.9213  P=0.615
#PERMANOVA_L_C Complexity: R2 = 0.01518 F=0.5705  P=0.671
#PERMANOVA_L_S Species: R2 = 0.50895 F = 3.8867 P = 0.001 ***


PERMANOVA2 = data.frame(PERMANOVA)
variable  = c('variable','variable','variable','variable')
PERMANOVA2$Significance <- "No Significant"
pval_thres <- 0.05
PERMANOVA2$Significance[which(PERMANOVA2$Pr..F. < pval_thres)] <- "Significant"
PERMANOVA2$Significance <- PERMANOVA2$Significance %>% factor

LC = data.frame(cbind(variable, row.names(PERMANOVA2), PERMANOVA2[,c(3,6)]))
colnames(LC) = c('variable','Effect', 'R2','Significance')

#Plot Variance
data_melt <- LC[-(4),]

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)#Turn your 'treatment' column into a character vector
data_melt$Effect <- factor(data_melt$Effect, levels=unique(rev(data_melt$Effect)))#Then turn it back into a factor with the levels in the correct order

mypal = c('white','red','blue','#528501','#2d4800','#1f4f58')

ggplot(data=data_melt, aes(x=variable, y=R2,fill=Effect)) +
  geom_bar(stat="identity",aes(color=Significance), size = 0.3,width = 0.9,linewidth=0.5) +# 
  scale_fill_manual(values = mypal) +
  theme_few() + #guides() + #color = FALSE, fill=FALSE
  labs(y="Variance explained (%)") +
  scale_color_manual(values = c('black',"red"),na.value =  "transparent")+
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0)),
        axis.text.x   = element_blank())

