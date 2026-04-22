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

#Control vs root exudates
ASV.50 <- names(sort(taxa_sums(ps.SR), TRUE)[1:50])
ps.SR.50 <- prune_taxa(ASV.50, ps.SR)

tableSyncom = cbind(ps.SR.50@sam_data, ps.SR.50@otu_table)
ncol(tableSyncom)
tableSyncom[1:7,1:7]
tableSyncom$Complexity = as.factor(tableSyncom$Complexity)
tableSyncom$Sample = as.factor(tableSyncom$Sample)

melted_tableSyncom <- tableSyncom %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

melted_tableSyncom_High = melted_tableSyncom[melted_tableSyncom$Root=='High',]
melted_tableSyncom_High_Control = melted_tableSyncom_High[melted_tableSyncom_High$Species=='Control',]
melted_tableSyncom_High_sControl = melted_tableSyncom_High[melted_tableSyncom_High$Species!='Control',]

melted_tableSyncom_Medium = melted_tableSyncom[melted_tableSyncom$Root=='Medium',]
melted_tableSyncom_Medium_Control = melted_tableSyncom_Medium[melted_tableSyncom_Medium$Species=='Control',]
melted_tableSyncom_Medium_sControl = melted_tableSyncom_Medium[melted_tableSyncom_Medium$Species!='Control',]

melted_tableSyncom_Low = melted_tableSyncom[melted_tableSyncom$Root=='Low',]
melted_tableSyncom_Low_Control = melted_tableSyncom_Low[melted_tableSyncom_Low$Species=='Control',]
melted_tableSyncom_Low_sControl = melted_tableSyncom_Low[melted_tableSyncom_Low$Species!='Control',]

log2foldchange_High <- c()
for(specie in melted_tableSyncom_High_sControl$Species  %>% unique){
  melted_sub2 <- melted_tableSyncom_High_sControl %>% subset(Species ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Control = data.frame(t(melted_tableSyncom_High_Control[,8]))
    High = data.frame(t(melted_sub[,8]))
    
    Control_mean = apply(Control, 1, mean) 
    High_mean = apply(High, 1, mean) 
    
    Control_High <- log2(High_mean+1) - log2(Control_mean+1) 
    
    Control_High_statistic = t.test(Control,High)
    pvalueControl_High = Control_High_statistic$p.value#Water_mean+1
    
    result = cbind(Control_High,pvalueControl_High)
    result = cbind(result,specie)
    row.names(result)= ASV
    log2foldchange_High <- rbind(log2foldchange_High,result)
  }
  Root <- data.frame(rep('High', 50))
  log2foldchange_High_Complete = cbind(Root,row.names(log2foldchange_High),log2foldchange_High)
  colnames(log2foldchange_High_Complete)=c('Root','ASV','diff','pvalue','Species_ID')
  
  log2foldchange_High_Complete <- rbind(log2foldchange_High_Complete,log2foldchange_High_Complete)
}

log2foldchange_Medium <- c()
for(specie in melted_tableSyncom_Medium_sControl$Species  %>% unique){
  melted_sub2 <- melted_tableSyncom_Medium_sControl %>% subset(Species ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Control = data.frame(t(melted_tableSyncom_High_Control[,8]))
    Medium = data.frame(t(melted_sub[,8]))
    
    Control_mean = apply(Control, 1, mean) 
    Medium_mean = apply(Medium, 1, mean) 
    
    High_Medium <- log2(Medium_mean+1) - log2(Control_mean+1) 
    
    High_Medium_statistic = t.test(Control,Medium)
    pvalueHigh_Medium = High_Medium_statistic$p.value#Water_mean+1
    
    result = cbind(High_Medium,pvalueHigh_Medium)
    result = cbind(result,specie)
    row.names(result)= ASV
    log2foldchange_Medium <- rbind(log2foldchange_Medium,result)
  }
  Root <- data.frame(rep('Medium', 50))
  log2foldchange_Medium_Complete = cbind(Root,row.names(log2foldchange_Medium),log2foldchange_Medium)
  colnames(log2foldchange_Medium_Complete)=c('Root','ASV','diff','pvalue','Species_ID')
  
  log2foldchange_Medium_Complete <- rbind(log2foldchange_Medium_Complete,log2foldchange_Medium_Complete)
}

log2foldchange_Low <- c()
for(specie in melted_tableSyncom_Low_sControl$Species  %>% unique){
  melted_sub2 <- melted_tableSyncom_Low_sControl %>% subset(Species ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Control = data.frame(t(melted_tableSyncom_High_Control[,8]))
    Low = data.frame(t(melted_sub[,8]))
    
    Control_mean = apply(Control, 1, mean) 
    Low_mean = apply(Low, 1, mean) 
    
    High_Low <- log2(Low_mean+1) - log2(Control_mean+1) 
    
    High_Low_statistic = t.test(Control,Low)
    pvalueHigh_Low = High_Low_statistic$p.value#Water_mean+1
    
    result = cbind(High_Low,pvalueHigh_Low)
    result = cbind(result,specie)
    row.names(result)= ASV
    log2foldchange_Low <- rbind(log2foldchange_Low,result)
  }
  Root <- data.frame(rep('Low', 50))
  log2foldchange_Low_Complete = cbind(Root,row.names(log2foldchange_Low),log2foldchange_Low)
  colnames(log2foldchange_Low_Complete)=c('Root','ASV','diff','pvalue','Species_ID')
  
  log2foldchange_Low_Complete <- rbind(log2foldchange_Low_Complete,log2foldchange_Low_Complete)
}


log2foldchange_Complete = rbind(log2foldchange_High_Complete,log2foldchange_Medium_Complete,log2foldchange_Low_Complete)

log2foldchange_Complete$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange_Complete$Significance[which(log2foldchange_Complete$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange_Complete$Significance <- log2foldchange_Complete$Significance %>% factor

log2foldchange_Complete[log2foldchange_Complete == '5.20554891117303',] <- 0
log2foldchange_Complete$diff [which(log2foldchange_Complete$diff  == '5.20554891117303')] <- "0"
log2foldchange_Complete$diff = as.numeric(log2foldchange_Complete$diff)

#Aggregate the dataframe to display as heatmap
display <- log2foldchange_Complete  %>%
  acast(formula = Species_ID*Root~ASV, mean,
        value.var = "diff")

display2=data.frame(display)
nrow(display2)

# Obtain the dendrogram
library(ggdendro)
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



log2foldchange_Complete$ASV = factor(log2foldchange_Complete$ASV, nutr_order)

log2foldchange_Complete$diff = as.numeric(log2foldchange_Complete$diff)
log2foldchange_Complete$Root = factor(log2foldchange_Complete$Root, c('Low','Medium','High'))
plt_hmap = ggplot(data = log2foldchange_Complete, aes(Species_ID,ASV)) + #
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.3,width = 0.9,height = 0.95) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Root, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-11,11),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs Control") + #Significance Genotype vs Col-0
  #theme(size_panel_border = 0.2)+
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05, size = 4),
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
