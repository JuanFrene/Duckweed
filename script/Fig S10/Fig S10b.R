packages <- c('ggthemes','dplyr', "ape", "Biostrings",  "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'Hmisc','reshape2')

sapply(packages, require, character.only = TRUE)              


setwd("G:/My Drive/labs/Nottingham/Duckweed/Metabolites 2 JGI/")
Metabolome_TP <- read.table("JGI targeted positive.txt", header = TRUE, row.names = 1)
Metabolome_TP2 = Metabolome_TP[Metabolome_TP$Exp!='Mono',]
Metabolome_TP2[,1:5]

melted_Metabolome_TP2 <- Metabolome_TP2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

head(melted_Metabolome_TP2)
melted_Metabolome_TP2$value = as.numeric(melted_Metabolome_TP2$value)

log2foldchange <- c()
for(specie in melted_Metabolome_TP2$Species  %>% unique){
  melted_sub <- melted_Metabolome_TP2 %>% subset(Species ==  specie) %>% droplevels
  for(var in melted_sub$variable  %>% unique){
    melted_sub2 <- melted_sub %>% subset(variable ==  var) %>% droplevels
    for(root in melted_sub2[melted_sub2$Treatment!='NB',]$Treatment  %>% unique){
      melted_sub3 <- melted_sub2[melted_sub2$Treatment!='NB',] %>% subset(Treatment ==  root) %>% droplevels #'ASV522'
      
      NB = data.frame(t(melted_sub2[melted_sub2$Treatment  == 'NB',][,5]))
      Compound = data.frame(t(melted_sub3[,5]))
      
      Control_mean = apply(NB, 1, mean) 
      Compound_mean = apply(Compound, 1, mean) 
      
      C_C_mean <- log2(Control_mean+1) - log2(Compound_mean+1) 
      
      C_C_mean_statistic = t.test(NB,Compound)
      pvalueB_C_C = C_C_mean_statistic$p.value
      
      result = cbind(C_C_mean,pvalueB_C_C)
      specie2 = specie
      row.names(result)= var
      result2 = cbind(specie2,result)
      result3 = cbind(root,result2)
      result4 = cbind(var,result3)
      log2foldchange <- rbind(log2foldchange,result4)
    }}}

log2foldchange2 = data.frame(log2foldchange)
head(log2foldchange2)
colnames(log2foldchange2) = c('Compound', 'Treatment', 'Species', 'diff', 'pvalue')

log2foldchange2$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange2$Significance[which(log2foldchange2$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange2$Significance <- log2foldchange2$Significance %>% factor
log2foldchange2$diff = as.numeric(log2foldchange2$diff)

head(log2foldchange2)
#Aggregate the dataframe to display as heatmap
display <- log2foldchange2  %>%
  acast(formula = Treatment~Compound, mean,
        value.var = "diff")

display2 = data.frame(t(display))

# Obtain the dendrogram
library(ggdendro)
dend <- as.dendrogram(hclust(dist(display2)))
dend_data <- dendro_data(dend)
ASV_order = dend_data$labels[,3]

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



log2foldchange2$Compound = factor(log2foldchange2$Compound, c(ASV_order))

log2foldchange2$diff[which(log2foldchange2$diff < -4)] <- -4
log2foldchange2$diff[which(log2foldchange2$diff > 4)] <- 4

plt_hmap = ggplot(data = log2foldchange2, aes(Treatment, Compound)) + 
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.9,height = 0.95) + #
  facet_grid(~Species, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-5,5),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs NB") + 
  #theme(size_panel_border = 0.2)7+
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05, size = 5),
        axis.text.y = element_text(size = 5),
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
plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.1, 1))

