packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')

sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")

# Control sequence ####
track <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/track.rds")
write.table(track, "track 16S.txt")

# 1. Metadata
meta <- read.table("Metadata 2.txt", header = TRUE, row.names = 1)

###ASV table
asv.table <- readRDS('G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/seqtab_final.rds')
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

##ASing to the syncom
taxaSyncom <- assignTaxonomy(asv.table2, "G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/Syncom_Sequence.pcr.unique.fasta", multithread=TRUE)#, minBoot = 98

#write.table(taxaSyncom, "taxaSyncom.txt")
#taxaSyncom = read.table("taxaSyncom.txt", header = TRUE, row.names = 1)

#Now we can make the phyloseq object
ps.Syncom <- phyloseq(asv.table2, tax_table(taxaSyncom), sample_data(meta))

taxa_names(ps.Syncom) <- paste0("ASV", seq(ntaxa(ps.Syncom)))

ps.Syncom = subset_taxa(ps.Syncom, Genus  != "NA")
#ps.SyncomSilva = subset_taxa(ps.SyncomSilva, Genus  != "NA")

ps.Syncom.3 = subset_samples(ps.Syncom, ID !=  "C-")
ps.pruned <- prune_taxa(taxa_sums(ps.Syncom.3)>=1, ps.Syncom.3)
ps.pruned3 <- prune_samples(sample_sums(ps.pruned)>1, ps.pruned) 
ps.perc <- transform_sample_counts(ps.pruned3, function(x) x / sum(x)) 


tableSyncom = cbind(ps.perc@sam_data, ps.perc@otu_table)
tableSyncom2=tableSyncom[,-c(1:6,10,11)][1:100]
tableSyncom3 = tableSyncom2[tableSyncom2$Drugs != 'DMSO',]

library(reshape2)
melted_table <- tableSyncom3 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

log2foldchange <- c()
for(specie in melted_table$Species  %>% unique){
  melted_sub <- melted_table %>% subset(Species ==  specie) %>% droplevels
  for(var in melted_sub$variable  %>% unique){
    melted_sub2 <- melted_sub %>% subset(variable ==  var) %>% droplevels #'ASV522'
    for(drug in melted_sub2[melted_sub2$Compound  !='Control',]$Drugs  %>% unique){
      melted_sub3 <- melted_sub2[melted_sub2$Compound  !='Control',] %>% subset(Drugs  ==  drug) %>% droplevels
      
      Control = data.frame(t(melted_sub2[melted_sub2$Compound  =='Control',][,5]))
      Compound = data.frame(t(melted_sub3[melted_sub3$Compound  =='Compound',][,5]))
      
      Control_mean = apply(Control, 1, mean) 
      Compound_mean = apply(Compound, 1, mean) 
      
      C_C_mean <- log2(Control_mean+1) - log2(Compound_mean+1) 
      
      C_C_mean_statistic = t.test(Control,Compound)
      pvalueB_C_C = C_C_mean_statistic$p.value
      C_C_mean_conf = t(data.frame(C_C_mean_statistic$conf.int))
      colnames(C_C_mean_conf) = c('inf','sup')
      
      result = cbind(C_C_mean,pvalueB_C_C)
      specie2 = specie
      row.names(result)= var
      result2 = cbind(specie2,result)
      result3 = cbind(drug,result2)
      result4 = cbind(var,result3)
      result5 = cbind(result4,C_C_mean_conf)
      log2foldchange <- rbind(log2foldchange,result5)
    }}}
log2foldchange2$Species
log2foldchange2 = data.frame(log2foldchange)
colnames(log2foldchange2) = c('ASV','Compound', 'Species', 'diff', 'pvalue', ' inf','sup')

log2foldchange2$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange2$Significance[which(log2foldchange2$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange2$Significance <- log2foldchange2$Significance %>% factor
log2foldchange2$diff = as.numeric(log2foldchange2$diff)

head(log2foldchange2)
#Aggregate the dataframe to display as heatmap
display <- log2foldchange2  %>%
  acast(formula = Species~ASV, mean,
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



log2foldchange2$ASV = factor(log2foldchange2$ASV, c(ASV_order))
log2foldchange2$Species = factor(log2foldchange2$Species, c('LY9250', "SP9505","SP7498","LP0049", "LJ9250","LP8539"))

plt_hmap = ggplot(data = log2foldchange2, aes(Compound, ASV)) + 
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.9,height = 0.95) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Species, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-0.25,0.25),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs Full") + #Significance Genotype vs Col-0
  #theme(size_panel_border = 0.2)+
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

plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.1, 1))

