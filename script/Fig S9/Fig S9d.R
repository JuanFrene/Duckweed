######Load packages
library(Rmisc)
library(reshape2)
library(ggrepel)
library(scales)
library(ggtree)
library(harrietr)
library(emmeans)
library(paletteer)
library(pals)
library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(ggdendro)



######Set work directory
setwd("G:/My Drive/labs/Nottingham/Duckweed/Drop out/Analysis files/")

######Open data file
Data <- read.table("LA7339 dropout.txt", header = TRUE, row.names = 1)
head(Data)
ID=rownames(Data)
Data= cbind(ID,Data)
Data$Syncom[which(Data$Syncom == 'Full_Syncom')] <- "Syncom"


### Now we need to summarize 
melted <- Data[,-(2:11)] %>% melt #Tab_f 
colnames(melted) <- c("ID","parameter","value")

###melted data and metadata
melted <- merge(Data[,c(1,4,10:11)],melted,  by='ID')

#Do a model per ion
Res_em <- NULL
Res_p <- NULL

for(batch in melted$Batch %>% unique){
  melted_batch <- melted %>% subset(Batch  == batch) %>% droplevels
  
  for(ion in melted_batch$variable %>% unique){
    melted_sub <- melted_batch %>% subset(variable  == ion) %>% droplevels
    m4 <- lm(data = melted_sub,
             formula = value ~ Syncom)
    m4_res <- emmeans(m4,pairwise~Syncom,ref =4, adjust = "none")
    m4_em <- m4_res$emmeans %>% as.data.frame
    m4_p <- m4_res$contrasts %>% as.data.frame
    m4_em$Ion <- ion
    m4_p$Ion <- ion
    Res_em <- rbind(Res_em,m4_em)
    Res_p <- rbind(Res_p,m4_p)
  }
}


#Arrange the pvalues dataframe
Res_p$contrast <- Res_p$contrast %>% gsub(pattern = " ",replacement = "")

#Now the contrast for Drop temperature
c9 <- Res_p$contrast %>% grep(pattern = "Syncom",value = T) %>% unique


####Subset the Drop temperature contrast
Res_p_Drop <- which(Res_p$contrast %in% c9) %>%
  Res_p[.,] %>% droplevels


###Correlation for the drop our differences.
#Aggregate the dataframe to display as heatmap
display_D <- Res_p_Drop  %>%
  acast(formula = Syncom~Ion, mean,
        value.var = "estimate")

display_D2 = data.frame(display_D)
head(display_D2)

##Corr y p-value
Corr020<-cor(display_D2[,-c(7:12)], method = 'pearson') #
#Corr020[is.na(Corr020)]<-0

p.mat020 <-cor_pmat(display_D2[,-c(7:12)], method = 'pearson') #[,-c(7:12)]
#p.mat020[is.na(p.mat020)]<-0


##Delimitando la matriz
CorrData5=Corr020
PvalueData5= p.mat020

sample_names <- colnames(CorrData5)

# Obtain the dendrogram
dend <- as.dendrogram(hclust(dist(CorrData5)))
dend_data <- dendro_data(dend)

dend_nutr <- as.dendrogram(hclust(dist(t(CorrData5))))
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
  mutate(x_center = (1:n()), 
         width = 1)

# Neglecting the gap parameters
heatmap_data <- CorrData5 %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

melted_CorrData5 <-CorrData5 %>% melt
melted_PvalueData5 <-PvalueData5 %>% melt
melted = cbind(melted_CorrData5, melted_PvalueData5[,3])
colnames(melted) = c('Var1','Var2','value', 'pvalue')
head(melted)

melted$Significance <- "No Significant"
pval_thres <- 0.05
melted$Significance[which(melted$pvalue < pval_thres)] <- "Significant"
melted$Significance <- melted$Significance %>% factor

melted$Var2 = factor(melted$Var2, nutr_order)

melted$Var1 = factor(melted$Var1, gene_pos_table$gene)


melterd2 = format(round(melted$value, 2), nsmall = 2) 
melted3 = cbind(melted,melterd2)

plt_hmap = ggplot(data = melted, aes(Var2,Var1)) +
  geom_raster(aes(fill = value))+
  theme_few() +
  geom_text(aes(label = melterd2), color = "black", size = 1.5) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.2,width = 0.9,height = 0.95) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-1,1),na.value = "#D9D9D9",name = "Corr") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significative") +
  labs(y='Microbial paramenters',x='Chemical parameters')+
  theme(axis.text.x = element_blank(), #element_text(angle = -45, hjust=-0.05
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
  theme(panel.grid.minor = element_blank(),axis.text.y=element_blank ())

plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.2, 1))