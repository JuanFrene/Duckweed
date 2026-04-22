###################################################################
library("FactoMineR")
library("factoextra")
library(ade4)
library(MASS)
library(ellipse)
library(ggplot2)
library(agricolae)
library(vegan)
library(ggthemes)
library(RVAideMemoire)
library(reshape2)
library(emmeans)
library(dplyr)
library(pals)
library(paletteer)

setwd("G:/My Drive/labs/Nottingham/Duckweed/exp Metabolites/Metabolome/")

Metabolome <- read.table("Table_metabolome.txt", header = TRUE, row.names = 1)
Metabolome1 = Metabolome[Metabolome$Treatment != 'QC',]
Metabolome2 = Metabolome1[Metabolome1$Species != 'LP8539',]

Metabolome2$Metabolite
ncol(Metabolome2)
Metabolome1[,1:5]

###log2foldchange
##transform our data into log2 base.
Metabolome_log = log2(Metabolome2[,-(1:5)])

##Add metadata
Metabolome1_log1 = (cbind(Metabolome2[,c(1:4)],Metabolome_log))
melted_Metabolome <- Metabolome1_log2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

list_log2foldchange <- c()

for(specie in melted_Metabolome$Species  %>% unique){
  melted_sub <- melted_Metabolome %>% subset(Species ==  specie) %>% droplevels
  
  for(id in melted_sub$variable      %>% unique){
    melted_sub2 <- melted_sub %>% subset(variable == id) %>% droplevels
    
    Syncom = data.frame(t(melted_sub2[melted_sub2$Treatment =='Syncom',][,6]))
    NB = data.frame(t(melted_sub2[melted_sub2$Treatment =='NB',][,6]))
    
    NB_mean = apply(NB, 1, mean) 
    Syncom_mean = apply(Syncom, 1, mean) 
    
    Syncommean_NBmean <- Syncom_mean - NB_mean 
    
    Syncommean_NBmean_statistic = t.test(Syncom,NB)
    pvalueSyncommean_NBmean = Syncommean_NBmean_statistic$p.value
    
    duckweed = data.frame(specie)
    Metabolite= data.frame(id)
    result = cbind(duckweed, Metabolite, Syncommean_NBmean,pvalueSyncommean_NBmean)
    
    
    list_log2foldchange <- rbind(list_log2foldchange,result)
  }}

head(list_log2foldchange)
list_log2foldchange2=data.frame(list_log2foldchange)

####Adjust the p-value
list_log2foldchange2$SignPvalue <- "NoSignificant"
pval_thres <- 0.05
list_log2foldchange2$SignPvalue[which(list_log2foldchange2$pvalueSyncommean_NBmean < pval_thres)] <- "q < 0.05"
list_log2foldchange2$SignPvalue <- list_log2foldchange2$SignPvalue %>% factor
colnames(list_log2foldchange2)=c('Species','Metabolite', 'log2foldchange', 'p.value', 'SignPval')
head(list_log2foldchange2)


###Ajustar por fdr
list_log2foldchange2$p.adj <- list_log2foldchange2$p.value %>% p.adjust(method = "fdr")
pval_thres2 <- 0.1
list_log2foldchange2$SignPajus <- "NoSignificant"
list_log2foldchange2$SignPajus[which(list_log2foldchange2$p.adj < pval_thres2)] <- "q < 0.1"
list_log2foldchange2$SignPajus <- list_log2foldchange2$SignPajus %>% factor
colnames(list_log2foldchange2)=c('Species','Metabolite', 'log2foldchange', 'p.value', 'SignPval','p.ajust','SignPajus')


###Clusterization
display <- list_log2foldchange2  %>%
  acast(formula = Species~Metabolite,fill = 0,
        value.var = "log2foldchange") %>%
  scale
display2=data.frame(display)

dend_nutr <- as.dendrogram(hclust(dist(t(display2))))
dend_nutr_data <- dendro_data(dend_nutr)
nutr_order = dend_nutr_data$labels[,3]

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_nutr_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_nutr_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = nutr_order) %>%
  mutate(x_center = (1:n()), width = 1)

# Neglecting the gap parameters
heatmap_data <- display2 %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  cross_join(gene_pos_table) %>%
  cross_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 1 * height))) + 0.1 * c(-1, 1) # extra spacing: 0.1

species_order = dend_nutr_data$labels[,3]
dend_nutr_data$labels
list_log2foldchange2$Metabolite=factor(list_log2foldchange2$Metabolite, nutr_order)
list_log2foldchange2$Species=factor(list_log2foldchange2$Species, c('SP7498','SI7820','LT9243','LA7339','SP9509','LY9205','SI9227','LJ9250','LP0049')) #species_order
nrow(list_log2foldchange2)

list_log2foldchange2$log2foldchange[which(list_log2foldchange2$log2foldchange > 0.05)] <- 0.05
list_log2foldchange2$log2foldchange[which(list_log2foldchange2$log2foldchange < -0.15)] <- -0.15

####Heatmap
plt_hmap = ggplot(data = list_log2foldchange2, aes(Species,Metabolite)) +
  geom_raster(aes(fill = log2foldchange))+
  theme_few() +
  geom_tile(aes(color = SignPval),fill = '#00000000', linewidth = 0.05,width = 0.9,height = 0.95) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-0.16,0.10),na.value = "#D9D9D9")+
  scale_color_manual(values = c('transparent',"black"),na.value =  "transparent",name = "Significance vs sterile medium") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05),axis.text.y=element_blank(),
        axis.title = element_blank())

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse() + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_tufte() + 
  theme(panel.grid.minor = element_blank(),axis.text.y=element_blank(),
        axis.title = element_blank()) #

library(cowplot)
plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.2, 1))

