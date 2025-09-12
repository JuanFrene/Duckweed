###PCoA Exudates
library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggplot2)
library(ggthemes)
library(ggdendro)
library(pals)
library(paletteer)
require(GGally)
require(CCA)
require(CCP)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exudados/")

Exudados <- read.table("Exudados.txt", header = TRUE, row.names = 1)
Exudados[,1:5]

Exudados$Subsample = as.factor(Exudados$Subsample)
Exudados$CCL = as.factor(Exudados$CCL)

melted_Exudados <- Exudados[1:250] %>% melt
gen.means <- aggregate(value~Species+variable, melted_Exudados, FUN=mean)
gen.means2 = cbind(gen.means,log10(gen.means[,3]))
colnames(gen.means2) = c('Species', 'variable', 'value','log10')

# Obtain the dendrogram
display <- gen.means2  %>%
  acast(formula = Species~variable, mean,
        value.var = "log10")

display2 = data.frame(display)


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



gen.means2$variable = factor(gen.means2$variable, c(nutr_order))#)
gen.means2$Species = factor(gen.means2$Species, dend_data$labels[,3])


plt_hmap = ggplot(data = gen.means2, aes(Species,variable)) +
  geom_raster(aes(fill = log10))+
  theme_few() +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",name = "Abundance") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs Water") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05),
        axis.text.y = element_text(size = 5))

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

