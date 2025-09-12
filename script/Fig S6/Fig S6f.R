###Heatmap
library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(ggdendro)
library(factoextra)
library(dplyr)
library(ggdendro)
library(reshape2)
library(cowplot)
library(ggthemes)
library(pals)
library(paletteer)


setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/Exp2/")

TablaCompleta = read.table("Root_trait_exp2_Complete.txt", header = TRUE, row.names = 1)
head(TablaCompleta)
TablaCompleta$Species = factor(TablaCompleta$Species, c('SP9509','SI7820','LJ9250','LP0049','LA7339','LP8539','LY5280')) 
TablaCompleta$ID = as.factor(TablaCompleta$ID)
TablaCompleta$Subsample = as.factor(TablaCompleta$Subsample)

melted_Table <- TablaCompleta[] %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

melted_Table_med = melted_Table[melted_Table$Nutrition!='Low',]
melted_Table_low = melted_Table[melted_Table$Nutrition!='Medium',]

log2foldchange_Table <- c()
#[melted_Table$variable!='cell_layer_n',]
for(treat in melted_Table$Microbiome  %>% unique){
  melted_sub <- melted_Table %>% subset(Microbiome  ==  treat) %>% droplevels
  for(specie in melted_sub$Species  %>% unique){
  melted_sub2 <- melted_sub %>% subset(Species  ==  specie) %>% droplevels
  for(var in melted_sub2$variable  %>% unique){
    melted_sub3 <- melted_sub2 %>% subset(variable ==  var) %>% droplevels
    for(nutrient in melted_sub3[melted_sub3$Nutrition !='High',]$Nutrition  %>% unique){
      melted_sub4 <- melted_sub3[melted_sub3$Nutrition !='High',] %>% subset(Nutrition  ==  nutrient) %>% droplevels
      
    
    NB = data.frame(t(melted_sub3[melted_sub3$Nutrition =='High',][,7]))
    B = data.frame(t(melted_sub4[,7]))
    
    B_mean = apply(B, 1, mean) 
    NB_mean = apply(NB, 1, mean) 
    
    B_NB_mean <- log2(B_mean+1) - log2(NB_mean+1)+0.1 
    
    B_NB_mean_statistic = t.test(B,NB)
    
    pvalueB_NB_mean = B_NB_mean_statistic$p.value
    
    result = cbind(B_NB_mean,pvalueB_NB_mean)
    specie2 = specie
    row.names(result)= var
    result2 = cbind(specie2,result)
    result3 = cbind(nutrient,result2)
    result4 = cbind(treat,result3)
    log2foldchange_Table <- rbind(log2foldchange_Table,result4)
  }}}}

nrow(log2foldchange_Table)

log2foldchange_Table_Complete = data.frame(cbind(row.names(log2foldchange_Table),log2foldchange_Table))
colnames(log2foldchange_Table_Complete)=c('Root_traits','Treatment','Nutrition','Species_ID','diff','pvalue')
log2foldchange_Table_Complete4 =log2foldchange_Table_Complete
head(log2foldchange_Table_Complete4)

log2foldchange_Table_Complete4$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange_Table_Complete4$Significance[which(log2foldchange_Table_Complete4$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange_Table_Complete4$Significance <- log2foldchange_Table_Complete4$Significance %>% factor

log2foldchange_Table_Complete4$var = paste(log2foldchange_Table_Complete4$Root_traits, log2foldchange_Table_Complete4$Nutrition, sep='_') 
log2foldchange_Table_Complete4$diff = as.numeric(log2foldchange_Table_Complete4$diff)

head(log2foldchange_Table_Complete4)

#Aggregate the dataframe to display as heatmap
display <- log2foldchange_Table_Complete4  %>%
  acast(formula = Species_ID~var, mean,
        value.var = "diff")

display2 = data.frame(display)


# Obtain the dendrogram
dend <- as.dendrogram(hclust(dist(display2)))
dend_data <- dendro_data(dend)

dend_nutr <- as.dendrogram(hclust(dist(t(display2))))
dend_nutr_data <- dendro_data(dend_nutr)
nutr_order = dend_nutr_data$labels[,3]

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = nutr_order) %>%
  group_by(x_center = (1:n()), width = 1)

# Neglecting the gap parameters
heatmap_data <- display2 %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  cross_join(gene_pos_table) %>%
  cross_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 1 * height))) + 0.1 * c(-1, 1) # extra spacing: 0.1



log2foldchange_Table_Complete4$var = factor(log2foldchange_Table_Complete4$var, c(nutr_order))
log2foldchange_Table_Complete4$Species_ID = factor(log2foldchange_Table_Complete4$Species_ID, dend_data$labels[,3])

log2foldchange_Table_Complete4$diff[which(log2foldchange_Table_Complete4$diff < -1)] <- "-1"
log2foldchange_Table_Complete4$diff[which(log2foldchange_Table_Complete4$diff > 1)] <- "1"

log2foldchange_Table_Complete4$diff = as.numeric(log2foldchange_Table_Complete4$diff)

plt_hmap = ggplot(data = log2foldchange_Table_Complete4, aes(var,Species_ID)) +
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.3,width = 0.9,height = 1) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Treatment, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-1,1),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('gray',"black"),na.value =  "transparent",name = "Significance vs Full nutrient") + #Significance Genotype vs Col-0
  #theme(size_panel_border = 0.2)+
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05, size = 5, color='black'),
        axis.title = element_blank(), legend.position = 'top')

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

plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.1, 1))


