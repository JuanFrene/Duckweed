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

Metabolome <- read.table("metabolite_log2foldchange2", header = TRUE, row.names = 1)

Species=c('SP7498','SI7820','LT9243','LA7339','SP9509','LY9205','SI9227','LJ9250','LP0049')
Complexity=c(7.25,6.8,6.2,5,7,6.8,7.6,5.33,5)
Met=data.frame(cbind(Species,Complexity))

TableMet6 = merge(list_log2foldchange2, Met, by= 'Species')


list_corrComplexity <- c()

overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

TableMet6$log2foldchange = as.numeric(TableMet6$log2foldchange)
TableMet6$Complexity = as.numeric(TableMet6$Complexity)

for(id in TableMet6$Metabolite  %>% unique){
  melted_sub <- TableMet6 %>% subset(Metabolite == id) %>% droplevels
  
  MetaCorr = lm(log2foldchange~Complexity, data=melted_sub)
  
  Metabolite= data.frame(id)
  Corr = cbind(Metabolite,overall_p(MetaCorr))
  
  list_corrComplexity <- rbind(list_corrComplexity,Corr)
}

list_corrComplexity2=data.frame(list_corrComplexity)
colnames(list_corrComplexity2)=c('Metabolite', 'p.value')

####Adjust the p-value
list_corrComplexity2$SignPvalue <- "NoSignificant"
pval_thres <- 0.05
list_corrComplexity2$SignPvalue[which(list_corrComplexity2$p.value < pval_thres)] <- "p < 0.05"
list_corrComplexity2$SignPvalue <- list_corrComplexity2$SignPvalue %>% factor

nrow(list_corrComplexity2)

###Subsample
Corr_metabolites <- list_corrComplexity2[list_corrComplexity2$SignPval ==  'p < 0.05',]
stelar = unique(Corr_metabolites$Metabolite)
Corr_metabolites = Corr_metabolites[with(Corr_metabolites, order(p.value)), ]


mms_all = merge(Corr_metabolites,list_log2foldchange2[,c(1:3,5)], by="Metabolite") 
head(mms_all)
colnames(mms_all) = c('Metabolite','p.value.corr','SignPvalue','Species','log2foldchange','p.value')
mms_all2 = merge(mms_all,Met, by="Species") 

nrow(mms_all)

###Clusterization
display_corr <- mms_all2  %>%
  acast(formula = Species~Metabolite,fill = 0,
        value.var = "log2foldchange") %>%
  scale


display_corr2=data.frame(display_corr)

dend_Species <- as.dendrogram(hclust(dist(display_corr2)))
dend_Species_data <- dendro_data(dend_Species)
Species_order = dend_Species_data$labels[,3]

dend_corr <- as.dendrogram(hclust(dist(t(display_corr2))))
dend_corr_data <- dendro_data(dend_corr)
corr_order = dend_corr_data$labels[,3]



# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data2 <- with(
  segment(dend_corr_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
corr_pos_table <- with(
  dend_corr_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_corr_table <- data.frame(sample = corr_order) %>%
  mutate(x_center = (1:n()), width = 1)

# Neglecting the gap parameters
heatmap_data2 <- display_corr2 %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  cross_join(corr_pos_table) %>%
  cross_join(sample_corr_table)

# Limits for the vertical axes
corr_axis_limits <- with(
  corr_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1) # extra spacing: 0.1

mms_all2$Metabolite=factor(mms_all2$Metabolite, corr_order)
mms_all2$Species=factor(mms_all2$Species, Species_order) #c('SP7498','SI7820','LT9243','LY9205','SP9509','LJ9250','LA7339','LP0049','SI9227')



####Heatmap
plt_hmap_corr = ggplot(data = mms_all2, aes(Species,Metabolite)) + #
  geom_raster(aes(fill = log2foldchange))+
  theme_few() +
  geom_tile(aes(color = p.value),fill = '#00000000', size = 0.2,width = 0.9,height = 0.95) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-0.16,0.05),na.value = "#D9D9D9")+
  scale_color_manual(values = c('transparent',"black"),na.value =  "transparent",name = "Significance vs Water") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05),axis.title = element_blank())

#Dendrogram plot
plt_dendr <- ggplot(segment_data2) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse() + 
  scale_y_continuous(breaks = corr_pos_table$y_center, 
                     labels = corr_pos_table$gene, 
                     limits = corr_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_tufte() + 
  theme(panel.grid.minor = element_blank(),axis.text.y=element_blank(),
        axis.title = element_blank()) #

plot_grid(plt_dendr, plt_hmap_corr, align = 'h', rel_widths = c(0.2, 1))

