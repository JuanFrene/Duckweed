packages <- c('ggthemes','dplyr', "ape", "Biostrings",  "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'Hmisc','reshape2', 'pals', 'paletteer')

sapply(packages, require, character.only = TRUE)              


setwd("G:/My Drive/labs/Nottingham/Duckweed/Metabolites 2 JGI/")

Metabolome_TN <- read.table("JGI targeted negative v2.txt", header = TRUE, row.names = 1)
#Metabolome_TN2 = Metabolome_TN[Metabolome_TN$Exp!='Mono',]
#head(Metabolome_TN2)
ncol(Metabolome_TN)
melted_Metabolome_TN2 <- Metabolome_TN[,-c(3,38,56,57)] %>% melt
melted_Metabolome_TN3 <- na.omit(melted_Metabolome_TN2)

Metabolome_TN_v <- melted_Metabolome_TN3  %>%
  acast(formula = Treatment*Species~variable, mean,
        value.var = "value")

#Metabolome_TN_v
zscoreN = scale(Metabolome_TN_v)
melted_zscoreN <- zscoreN[,-c(35,54,53)] %>% melt

Metabolome_TN3 = cbind(melted_zscoreN, melted_zscoreN[,1])
colnames(Metabolome_TN3) = c('Species', 'parameter', 'value','Treatment')
head(Metabolome_TN3)

Species <- Metabolome_TN3$Species %>% gsub(pattern = "Syncom_",replacement = "") %>%
  gsub(pattern = "NB_", replacement = "") %>%
  gsub(pattern = "dA245_", replacement = "") %>%
  gsub(pattern = "dA289_", replacement = "") %>%
  gsub(pattern = "dP65_", replacement = "") %>%
  gsub(pattern = "mA245_", replacement = "") %>%
  gsub(pattern = "mA289_", replacement = "") %>%
  gsub(pattern = "mP65_", replacement = "") %>%
  as.data.frame

head(Metabolome_TN3)

Treatment <- Metabolome_TN3$Treatment %>% gsub(pattern = "_LA7339",replacement = "") %>%
  gsub(pattern = "_SP9509", replacement = "") %>%
  as.data.frame

Metabolome_TN4= cbind(Species, Treatment, Metabolome_TN3[,2:3])
colnames(Metabolome_TN4) = c('Species','Treatment', 'variable', 'value')
head(Metabolome_TN4)

#Metabolome_TN5 = Metabolome_TN4[Metabolome_TN4$Treatment != 'mP65',] 
#Metabolome_TN6 = Metabolome_TN5[Metabolome_TN5$Treatment != 'mA245',] 
#Metabolome_TN7 = Metabolome_TN6[Metabolome_TN6$Treatment != 'mA289',] 
#Do a model per ion
head(melted_Metabolome_TN3)
TN_em <- NULL
TN_p <- NULL

#for(species in melted_Metabolome_TN3$Species %>% unique){
melted_Metabolome_TN_SP = melted_Metabolome_TN3[melted_Metabolome_TN3$Species== 'SP9509',]
for(subs in melted_Metabolome_TN_SP$variable %>% unique){
  melted_sub <- melted_Metabolome_TN_SP %>% subset(variable  == subs) %>% droplevels
  m4 <- lm(data = melted_sub,
           formula = value ~ Treatment)
  m4_res <- emmeans(m4,pairwise~Treatment,ref =4, adjust = "none")
  m4_em <- m4_res$emmeans %>% as.data.frame
  m4_p <- m4_res$contrasts %>% as.data.frame
  m4_em$Species <- 'SP9509'
  m4_em$Metabolite <- subs
  m4_p$Species <- 'SP9509'
  m4_p$Metabolite <- subs
  TN_em <- rbind(TN_em,m4_em)
  TN_p <- rbind(TN_p,m4_p)
}
SP_N = TN_p

melted_Metabolome_TN_LA = melted_Metabolome_TN3[melted_Metabolome_TN3$Species== 'LA7339',]
for(subs in melted_Metabolome_TN_LA$variable %>% unique){
  melted_sub <- melted_Metabolome_TN_LA %>% subset(variable  == subs) %>% droplevels
  m4 <- lm(data = melted_sub,
           formula = value ~ Treatment)
  m4_res <- emmeans(m4,pairwise~Treatment,ref =4, adjust = "none")
  m4_em <- m4_res$emmeans %>% as.data.frame
  m4_p <- m4_res$contrasts %>% as.data.frame
  m4_em$Species <- 'LA7339'
  m4_em$Metabolite <- subs
  m4_p$Species <- 'LA7339'
  m4_p$Metabolite <- subs
  TN_em <- rbind(TN_em,m4_em)
  TN_p <- rbind(TN_p,m4_p)
}

TN_p2 =rbind(TN_p,SP_N)

unique(TN_p2$Species)

#Arrange the pvalues dataframe
TN_p2$contrast <- TN_p2$contrast %>% gsub(pattern = " ",replacement = "")

#Now the contrast for Room temperature
c9 <- TN_p2$contrast %>% grep(pattern = "NB",value = T) %>% unique

####Subset the contrast
TN_p_com <- which(TN_p2$contrast %in% c9) %>%
  TN_p2[.,] %>% droplevels

nrow(TN_p_com)

TN_p_com$Significance <- "NoSignificant"
pval_thres <- 0.1
TN_p_com$Significance[which(TN_p_com$p.value < pval_thres)] <- "q < 0.1"
TN_p_com$Significance <-TN_p_com$Significance %>% factor

TN_p_com$contrast <-TN_p_com$contrast %>% gsub(pattern = "-NB",replacement = "")
TN_p_com$contrast <-TN_p_com$contrast %>% gsub(pattern = "NB-",replacement = "")

TN_p_com$Id <- paste(TN_p_com$Species,TN_p_com$Metabolite,TN_p_com$contrast,sep = "_")
Metabolome_TN4$Id <- paste(Metabolome_TN4$Species,Metabolome_TN4$variable,Metabolome_TN4$Treatment,sep = "_")

Metabolome_TN5 = merge(Metabolome_TN4, TN_p_com[,9:10], by = 'Id')

TN_p_com$Significance <- "NoSignificant"
pval_thres <- 0.1
TN_p_com$Significance[which(TN_p_com$p.value < pval_thres)] <- "q < 0.1"
TN_p_com$Significance <-TN_p_com$Significance %>% factor

head(Metabolome_TN5)
head(Metabolome_TN4)

Metabolome_TN4$Significance <- "NoSignificant"

Metabolome_TN6 = rbind(Metabolome_TN5,
                       Metabolome_TN4[Metabolome_TN4$Treatment=='NB',],
                       Metabolome_TN4[Metabolome_TN4$variable=='X4.Coumaric_acid',])


#Aggregate the dataframe to display as heatmap
display <- Metabolome_TN4  %>%
  acast(formula = Treatment~variable, mean,
        value.var = "value")

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



Metabolome_TN6$variable = factor(Metabolome_TN6$variable, c(ASV_order))
Metabolome_TN6$Treatment = factor(Metabolome_TN6$Treatment, c('NB',"Syncom","dA289","dA245","dP65","mA289","mA245","mP65"))

plt_hmap = ggplot(data = Metabolome_TN6, aes(Treatment, variable)) + 
  geom_raster(aes(fill = value)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.9,height = 0.95) + #
  facet_grid(~Species, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-3,3),na.value = "black",name = "Abundance") +
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
Negative = plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.2, 1))


Metabolome_TP <- read.table("JGI targeted positive v2.txt", header = TRUE, row.names = 1)
head(Metabolome_TP)
melted_Metabolome_TP2 <- Metabolome_TP[,-c(3,28,31)] %>% melt
melted_Metabolome_TP3 <- na.omit(melted_Metabolome_TP2)

Metabolome_TP_v <- melted_Metabolome_TP3  %>%
  acast(formula = Treatment*Species~variable, mean,
        value.var = "value")

#Metabolome_TP_v
zscoreN = scale(Metabolome_TP_v)
melted_zscoreN <- zscoreN %>% melt

Metabolome_TP3 = cbind(melted_zscoreN, melted_zscoreN[,1])
colnames(Metabolome_TP3) = c('Species', 'parameter', 'value','Treatment')
head(Metabolome_TP3)

Species <- Metabolome_TP3$Species %>% gsub(pattern = "Syncom_",replacement = "") %>%
  gsub(pattern = "NB_", replacement = "") %>%
  gsub(pattern = "dA245_", replacement = "") %>%
  gsub(pattern = "dA289_", replacement = "") %>%
  gsub(pattern = "dP65_", replacement = "") %>%
  gsub(pattern = "mA245_", replacement = "") %>%
  gsub(pattern = "mA289_", replacement = "") %>%
  gsub(pattern = "mP65_", replacement = "") %>%
  as.data.frame

head(Metabolome_TP3)

Treatment <- Metabolome_TP3$Treatment %>% gsub(pattern = "_LA7339",replacement = "") %>%
  gsub(pattern = "_SP9509", replacement = "") %>%
  as.data.frame

Metabolome_TP4= cbind(Species, Treatment, Metabolome_TP3[,2:3])
colnames(Metabolome_TP4) = c('Species','Treatment', 'variable', 'value')
Metabolome_TP4$Significance <- "NoSignificant"
head(Metabolome_TP4)

#Do a model per ion
head(melted_Metabolome_TP3)
TP_em <- NULL
TP_p <- NULL

#for(species in melted_Metabolome_TN3$Species %>% unique){
melted_Metabolome_TP_SP = melted_Metabolome_TP3[melted_Metabolome_TP3$Species== 'SP9509',]
for(subs in melted_Metabolome_TP_SP$variable %>% unique){
  melted_sub <- melted_Metabolome_TP_SP %>% subset(variable  == subs) %>% droplevels
  m4 <- lm(data = melted_sub,
           formula = value ~ Treatment)
  m4_res <- emmeans(m4,pairwise~Treatment,ref =4, adjust = "none")
  m4_em <- m4_res$emmeans %>% as.data.frame
  m4_p <- m4_res$contrasts %>% as.data.frame
  m4_em$Species <- 'SP9509'
  m4_em$Metabolite <- subs
  m4_p$Species <- 'SP9509'
  m4_p$Metabolite <- subs
  TP_em <- rbind(TP_em,m4_em)
  TP_p <- rbind(TP_p,m4_p)
}

melted_Metabolome_TP_LA = melted_Metabolome_TP3[melted_Metabolome_TP3$Species == 'LA7339',]
for(subs in melted_Metabolome_TP_LA$variable %>% unique){
  melted_sub <- melted_Metabolome_TP_LA %>% subset(variable  == subs) %>% droplevels
  m4 <- lm(data = melted_sub,
           formula = value ~ Treatment)
  m4_res <- emmeans(m4,pairwise~Treatment,ref =4, adjust = "none")
  m4_em <- m4_res$emmeans %>% as.data.frame
  m4_p <- m4_res$contrasts %>% as.data.frame
  m4_em$Species <- 'LA7339'
  m4_em$Metabolite <- subs
  m4_p$Species <- 'LA7339'
  m4_p$Metabolite <- subs
  TP_em <- rbind(TP_em,m4_em)
  TP_p <- rbind(TP_p,m4_p)
}


unique(TP_p$Species)

#Arrange the pvalues dataframe
TP_p$contrast <- TP_p$contrast %>% gsub(pattern = " ",replacement = "")

#Now the contrast for Room temperature
c9 <- TP_p$contrast %>% grep(pattern = "NB",value = T) %>% unique

####Subset the contrast
TP_p_com <- which(TP_p$contrast %in% c9) %>%
  TP_p[.,] %>% droplevels

nrow(TP_p_com)

TP_p_com$Significance <- "NoSignificant"
pval_thres <- 0.1
TP_p_com$Significance[which(TP_p_com$p.value < pval_thres)] <- "q < 0.1"
TP_p_com$Significance <-TP_p_com$Significance %>% factor

TP_p_com$contrast <-TP_p_com$contrast %>% gsub(pattern = "-NB",replacement = "")
TP_p_com$contrast <-TP_p_com$contrast %>% gsub(pattern = "NB-",replacement = "")

TP_p_com$Id <- paste(TP_p_com$Species,TP_p_com$Metabolite,TP_p_com$contrast,sep = "_")
Metabolome_TP4$Id <- paste(Metabolome_TP4$Species,Metabolome_TP4$variable,Metabolome_TP4$Treatment,sep = "_")

Metabolome_TP5 = merge(Metabolome_TP4[,-(5)], TP_p_com[,9:10], by = 'Id')

head(Metabolome_TP5)
head(Metabolome_TP4)

Metabolome_TP6 = rbind(Metabolome_TP5,
                       Metabolome_TP4[Metabolome_TP4$Treatment=='NB',])


#Aggregate the dataframe to display as heatmap
display <- Metabolome_TP4  %>%
  acast(formula = Treatment~variable, mean,
        value.var = "value")

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

Metabolome_TP6$variable = factor(Metabolome_TP6$variable, c(ASV_order))
Metabolome_TP6$Treatment = factor(Metabolome_TP6$Treatment, c('NB',"Syncom","dA289","dA245","dP65","mA289","mA245","mP65"))

plt_hmap = ggplot(data = Metabolome_TP6, aes(Treatment, variable)) + 
  geom_raster(aes(fill = value)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.9,height = 0.95) + #
  facet_grid(~Species, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-3,3),na.value = "black",name = "Abundance") +
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

#library(cowplot)
Positive = plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.2, 1))

plot_grid(Negative, Positive, align = 'h', rel_widths = c(0.46, 0.54),   labels = "auto")

