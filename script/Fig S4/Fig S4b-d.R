packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest",'dada2', "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", "reshape2")
sapply(packages, require, character.only = TRUE)              

library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(harrietr)
library(egg)
library(paletteer)
library(ggtree)
library(scales)
library(car)
library(Rmisc)
library(factoextra)
library(ggdendro)
library(cowplot)
library(ggcorrplot)


setwd("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/")

###Track sequencing abundance
#track <- readRDS("C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/track240-200.rds")
#write.table(track, "track 16S 240-240 Exp1.txt")

# 1. Taxonomy Table
taxa <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/tax_final250_220.rds")
# check your taxonomic classifications #
taxa.print<- taxa
rownames(taxa.print)
head(taxa.print)
unname(head(taxa))

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)
#write.table(ps2@sam_data, file = "C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/metadataP.tsv", row.names=FALSE, sep="\t")

###ASV table
asv.table<- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/seqtab_final250_220.rds")
#write.table(asv.table, file = "C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/asv_table.txt", row.names=FALSE, sep="\t")
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)


# 14. Assign taxonomy gg_13_8 ####
#otumat <- read.table("clipboard", header = TRUE)
#taxa.print2<- otumat
#rownames(taxa.print2)
#head(taxa.print2)
#unname(head(otumat))

# check your taxonomic classifications #
#taxgg <- assignTaxonomy(asv.table2, "C:/Users/juanp/Downloads/gg_13_8_train_set_97.fa.gz", multithread=TRUE)

#Now we can make the phyloseq object
ps2 <- phyloseq(asv.table2, tax_table(taxa), sample_data(meta))
#ps2gg <- phyloseq(asv.table2, tax_table(taxgg), sample_data(meta))

#Random phylogenetic tree
random_tree = rtree(ntaxa(ps2), rooted=TRUE, tip.label=taxa_names(ps2))

ps.1 = merge_phyloseq(ps2, random_tree)

# while it doesn't seem to be the case, if it were we would use the following lines to remove samples with low sequence counts
set.seed(500)
ps.2 = subset_taxa(ps.1, Kingdom == "Bacteria" | Kingdom == "Archaea")
ps.pruned <- prune_samples(sample_sums(ps.2)>=2, ps.2)
ps.pruned_taxa <- filter_taxa(ps.2, function(x) sum(x) > .005, TRUE)

ps.3 = subset_samples(ps.pruned, Compartment !=  "BLANK")
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))

ps.perc.3 <- transform_sample_counts(ps.3, function(x) x / sum(x)) 


ps.3.Root = subset_samples(ps.3, Compartment ==  "Root")
ps.3.Water = subset_samples(ps.3, Compartment ==  "Water")
ps.3.SInoculum = subset_samples(ps.3, Compartment !=  "Inoculum")
ps.3.SWater = subset_samples(ps.3.SInoculum, Compartment !=  "Water")
ps.3.Front = subset_samples(ps.3, Compartment ==  "Front")
ps.3.SC = subset_samples(ps.3, Specie !=  "Control")
ps.3.SCSI = subset_samples(ps.3.SC, Specie !=  "Inoculum")

ps.3.SWater.perc<- transform_sample_counts(ps.3.SWater, function(x) x / sum(x))
ps.perc.Root = subset_samples(ps.3.SWater.perc, Compartment ==  "Root")
ps.perc.Front = subset_samples(ps.3.SWater.perc, Compartment ==  "Front")


#Root and Front vs water
ps.3.SInoculum
ps.3.SInoculum_core = core(ps.3.SInoculum , detection = 0.001, prevalence = 20/100)
ps.Syncom.Syncom31 = subset_samples(ps.3.SInoculum.H_core, Compartment !=  "Root")
ps.Syncom.Syncom32 = subset_samples(ps.3.SInoculum.H_core, Specie !=  "Control")

tableSyncom = cbind(ps.Syncom.Syncom32@sam_data, ps.Syncom.Syncom32@otu_table)
tableSyncom2=tableSyncom[,-c(1:7,10,13)]
tableSyncom2$Complexity = as.factor(tableSyncom2$Complexity)
tableSyncom2$Subsample = as.factor(tableSyncom2$Subsample)

melted_tableSyncom <- tableSyncom2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

melted_tableSyncom_High

log2foldchange_Front <- c()
log2foldchange_Front_Complete2 =c()

for(specie in melted_tableSyncom$Specie.1  %>% unique){
  melted_sub2 <- melted_tableSyncom %>% subset(Specie.1 ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Front = data.frame(t(melted_sub[melted_sub$Compartment=='Front',][,8]))
    Water = data.frame(t(melted_sub[melted_sub$Compartment=='Water',][,8]))
    
    Front_mean = apply(Front, 1, mean) 
    Water_mean = apply(Water, 1, mean) 
    
    Frontmean_Watermean <- log2(Front_mean+1) - log2(Water_mean+1) 
    
    Frontmean_Watermean_statistic = t.test(Front,Water)
    pvalueFrontmean_Watermean = Frontmean_Watermean_statistic$p.value#Water_mean+1
    
    result = cbind(Frontmean_Watermean,pvalueFrontmean_Watermean)
    row.names(result)= ASV
    log2foldchange_Front <- rbind(log2foldchange_Front,result)
  }
  Compartment <- data.frame(rep('Front', 95))
  Species_ID= data.frame(rep(specie, 95))
  log2foldchange_Front_Complete = cbind(Compartment,Species_ID,row.names(log2foldchange_Front),log2foldchange_Front)
  colnames(log2foldchange_Front_Complete)=c('Compartment','Species_ID','ASV','diff','pvalue')
  
  log2foldchange_Front_Complete2 <- rbind(log2foldchange_Front_Complete2,log2foldchange_Front_Complete)
}

log2foldchange_Root <- c()
log2foldchange_Root_Complete2 =c()

for(specie in melted_tableSyncom$Specie.1  %>% unique){
  melted_sub2 <- melted_tableSyncom %>% subset(Specie.1 ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Root = data.frame(t(melted_sub[melted_sub$Compartment=='Root',][,7]))
    Water = data.frame(t(melted_sub[melted_sub$Compartment=='Water',][,7]))
    
    Root_mean = apply(Root, 1, mean) 
    Water_mean = apply(Water, 1, mean) 
    
    Rootmean_Watermean <- log2(Root_mean+1) - log2(Water_mean+1) 
    
    Rootmean_Watermean_statistic = t.test(Root,Water)
    pvalueRootmean_Watermean = Rootmean_Watermean_statistic$p.value#Water_mean+1
    
    result = cbind(Rootmean_Watermean,pvalueRootmean_Watermean)
    row.names(result)= ASV
    log2foldchange_Root <- rbind(log2foldchange_Root,result)
  }
  Compartment <- data.frame(rep('Root', 97))
  Species_ID= data.frame(rep(specie, 97))
  log2foldchange_Root_Complete = cbind(Compartment,Species_ID,row.names(log2foldchange_Root),log2foldchange_Root)
  colnames(log2foldchange_Root_Complete)=c('Compartment','Species_ID','ASV','diff','pvalue')
  
  log2foldchange_Root_Complete2 <- rbind(log2foldchange_Root_Complete2,log2foldchange_Root_Complete)
}

log2foldchange_Complete = rbind(log2foldchange_Root_Complete2,log2foldchange_Front_Complete2)

log2foldchange_Complete$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange_Complete$Significance[which(log2foldchange_Complete$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange_Complete$Significance <- log2foldchange_Complete$Significance %>% factor

log2foldchange_Complete$Compartment [which(log2foldchange_Complete$Compartment  == 'Front')] <- "Frond"


#Aggregate the dataframe to display as heatmap

display <- log2foldchange_Complete  %>%
  acast(formula = Species_ID*Compartment~ASV,fill = 0,
        value.var = "diff") %>% scale

display2=data.frame(display)

library(factoextra)
res.display2 <- hclust(dist(t(tableSyncom2[,-(1:5)])),  method = "ward.D")#[c(1,7,13,19,25,31,37),]
fviz_dend(res.display2, cex = 0.35, k = 3, palette = "jco",rect = TRUE, main="Cluster by Samples",  horiz = TRUE) 

log2foldchange_Complete$ASV=factor(log2foldchange_Complete$ASV, label_order_vswater)
log2foldchange_Complete$Compartment=factor(log2foldchange_Complete$Compartment, c('Root','Frond'))

ggplot(data = log2foldchange_Complete, aes(Species_ID,ASV)) +
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.2,width = 1,height = 1) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Compartment, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-0.4,0.4),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs Water") + #Significance Genotype vs Col-0
  #theme(size_panel_border = 0.2)+
  theme(axis.text.x = element_text(angle = -90,vjust=-0.05,hjust=-0.05, size = 7)) #, 



#Root and Front vs water

ps.phyla.perc <- taxa_level(ps.3.SInoculum, "Genus")

ps.3.SInoculum.H_core_phylum = core(ps.phyla.perc, detection = 0.00001, prevalence = 20/100)
ps.Syncom.Syncom31ps.phyla.perc = subset_samples(ps.3.SInoculum.H_core_phylum, Compartment !=  "Root")
ps.Syncom.Syncom32 = subset_samples(ps.3.SInoculum.H_core_phylum, Specie !=  "Control")

tableSyncom = cbind(ps.Syncom.Syncom32@sam_data, ps.Syncom.Syncom32@otu_table)
tableSyncom2=tableSyncom[,-c(1:8,10,14,15)]
tableSyncom2$Complexity = as.factor(tableSyncom2$Complexity)
tableSyncom2$Subsample = as.factor(tableSyncom2$Subsample)

melted_tableSyncom <- tableSyncom2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

log2foldchange_Front <- c()
log2foldchange_Front_Complete2 =c()

for(specie in melted_tableSyncom$Specie.1  %>% unique){
  melted_sub2 <- melted_tableSyncom %>% subset(Specie.1 ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Front = data.frame(t(melted_sub[melted_sub$Compartment=='Front',][,6]))
    Water = data.frame(t(melted_sub[melted_sub$Compartment=='Water',][,6]))
    
    Front_mean = apply(Front, 1, mean) 
    Water_mean = apply(Water, 1, mean) 
    
    Frontmean_Watermean <- log2(Front_mean+1) - log2(Water_mean+1) 
    
    Frontmean_Watermean_statistic = t.test(Front,Water)
    pvalueFrontmean_Watermean = Frontmean_Watermean_statistic$p.value#Water_mean+1
    
    result = cbind(Frontmean_Watermean,pvalueFrontmean_Watermean)
    row.names(result)= ASV
    log2foldchange_Front <- rbind(log2foldchange_Front,result)
  }
  Compartment <- data.frame(rep('Front', 53))
  Species_ID= data.frame(rep(specie, 53))
  log2foldchange_Front_Complete = cbind(Compartment,Species_ID,row.names(log2foldchange_Front),log2foldchange_Front)
  colnames(log2foldchange_Front_Complete)=c('Compartment','Species_ID','ASV','diff','pvalue')
  
  log2foldchange_Front_Complete2 <- rbind(log2foldchange_Front_Complete2,log2foldchange_Front_Complete)
}

log2foldchange_Root <- c()
log2foldchange_Root_Complete2 =c()

for(specie in melted_tableSyncom$Specie.1  %>% unique){
  melted_sub2 <- melted_tableSyncom %>% subset(Specie.1 ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Root = data.frame(t(melted_sub[melted_sub$Compartment=='Root',][,6]))
    Water = data.frame(t(melted_sub[melted_sub$Compartment=='Water',][,6]))
    
    Root_mean = apply(Root, 1, mean) 
    Water_mean = apply(Water, 1, mean) 
    
    Rootmean_Watermean <- log2(Root_mean+1) - log2(Water_mean+1) 
    
    Rootmean_Watermean_statistic = t.test(Root,Water)
    pvalueRootmean_Watermean = Rootmean_Watermean_statistic$p.value#Water_mean+1
    
    result = cbind(Rootmean_Watermean,pvalueRootmean_Watermean)
    row.names(result)= ASV
    log2foldchange_Root <- rbind(log2foldchange_Root,result)
  }
  Compartment <- data.frame(rep('Root', 53))
  Species_ID= data.frame(rep(specie, 53))
  log2foldchange_Root_Complete = cbind(Compartment,Species_ID,row.names(log2foldchange_Root),log2foldchange_Root)
  colnames(log2foldchange_Root_Complete)=c('Compartment','Species_ID','ASV','diff','pvalue')
  
  log2foldchange_Root_Complete2 <- rbind(log2foldchange_Root_Complete2,log2foldchange_Root_Complete)
}

log2foldchange_Complete = rbind(log2foldchange_Root_Complete2,log2foldchange_Front_Complete2)

log2foldchange_Complete$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange_Complete$Significance[which(log2foldchange_Complete$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange_Complete$Significance <- log2foldchange_Complete$Significance %>% factor
log2foldchange_Complete$Compartment [which(log2foldchange_Complete$Compartment  == 'Front')] <- "Frond"
#log2foldchange_Complete$ASV [which(log2foldchange_Complete$ASV  == 'LiUU-11-161')] <- "LiUU_11_161"
#log2foldchange_Complete$ASV [which(log2foldchange_Complete$ASV  == 'NS11-12_marine_group')] <- "NS11_12_marine_group"
log2foldchange_Complete$ASV [which(log2foldchange_Complete$ASV  == 'TM6_(Dependentiae)')] <- "TM6"

#Aggregate the dataframe to display as heatmap

display <- log2foldchange_Complete[log2foldchange_Complete$Compartment=='Root',]  %>%
  acast(formula = Species_ID~ASV,mean,
        value.var = "diff")

display2=data.frame(display)


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



log2foldchange_Complete$ASV = factor(log2foldchange_Complete$ASV, nutr_order)
log2foldchange_Complete$Compartment = factor(log2foldchange_Complete$Compartment, c('Root', 'Frond'))

log2foldchange_Complete$Species_ID = factor(log2foldchange_Complete$Species_ID, dend_data$labels[,3])

plt_hmap = ggplot(data = log2foldchange_Complete, aes(ASV,Species_ID)) +
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.9,height = 0.95) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Compartment, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-11,11),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs Water") + #Significance Genotype vs Col-0
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

plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.1, 1))
