packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')
sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/")

# Control sequence ####
track240_240 <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/track240240.rds")

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

###ASV table
asv.table <- readRDS('G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/seqtab_final240_220.rds')
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

##ASing to the syncom
taxaSyncom <- assignTaxonomy(asv.table2, "G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/Syncom_Sequence.pcr.unique.fasta", multithread=TRUE)#, minBoot = 98

#Now we can make the phyloseq object
ps.Syncom <- phyloseq(asv.table2, tax_table(taxaSyncom), sample_data(meta))
taxa_names(ps.Syncom) <- paste0("ASV", seq(ntaxa(ps.Syncom)))

ps.Syncom = subset_taxa(ps.Syncom, Genus  != "NA")

ps.Syncom.3 = subset_samples(ps.Syncom, Treatment !=  "ConNeg")
ps.pruned <- prune_taxa(taxa_sums(ps.Syncom.3)>=1, ps.Syncom.3)

ps.Syncom.3.H <- transform(ps.pruned, "hellinger")
ps.Syncom.3.perc <- transform_sample_counts(ps.pruned, function(x) x / sum(x)) 

ps.Syncom.SC = subset_samples(ps.pruned, Species !=  "Control")
ps.Syncom.SC2 = subset_samples(ps.Syncom.SC, Species !=  "Inoculum")

ps.Syncom.NOSyncom = subset_samples(ps.Syncom.SC2, Treatment !=  "SynCom")
ps.Syncom.NOSyncom2 <- prune_taxa(taxa_sums(ps.Syncom.NOSyncom)>0, ps.Syncom.NOSyncom)

ps.Syncom.Syncom = subset_samples(ps.Syncom.SC2, Treatment ==  "SynCom")
ps.Syncom.Syncom2 <- prune_taxa(taxa_sums(ps.Syncom.Syncom)>1, ps.Syncom.Syncom)

ps.Syncom.SyncomW = subset_samples(ps.Syncom.Syncom, Compartment !=  "Water")
ps.Syncom.SyncomW2 <- prune_taxa(taxa_sums(ps.Syncom.SyncomW)>440, ps.Syncom.SyncomW)

tableSyncom = cbind(ps.Syncom.SyncomW2@sam_data, ps.Syncom.SyncomW2@otu_table)
tableSyncom2=tableSyncom[,-c(1,2,3,8)]
tableSyncom2$cell_layer_n = as.factor(tableSyncom2$cell_layer_n)

melted_tableSyncom <- tableSyncom2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

melted_tableSyncom_med = melted_tableSyncom[melted_tableSyncom$Nutrient!='Low',]
melted_tableSyncom_Low = melted_tableSyncom[melted_tableSyncom$Nutrient!='Medium',]

log2foldchange_Medium <- c()
for(part in melted_tableSyncom_med$Compartment  %>% unique){
  melted_sub <- melted_tableSyncom_med %>% subset(Compartment  ==  part) %>% droplevels
  for(specie in melted_sub$Species  %>% unique){
    melted_sub2 <- melted_sub %>% subset(Species  ==  specie) %>% droplevels
    for(var in melted_sub2$variable  %>% unique){
      melted_sub3 <- melted_sub2 %>% subset(variable ==  var) %>% droplevels #'ASV522'
      
      NB = data.frame(t(melted_sub3[melted_sub3$Nutrient  =='High',][,7]))
      B = data.frame(t(melted_sub3[melted_sub3$Nutrient  =='Medium',][,7]))
      
      B_mean = apply(B, 1, mean) 
      NB_mean = apply(NB, 1, mean) 
      
      B_NB_mean <- log2(B_mean+1) - log2(NB_mean+1) 
      
      B_NB_mean_statistic = t.test(B,NB)
      
      pvalueB_NB_mean = B_NB_mean_statistic$p.value
      
      result = cbind(B_NB_mean,pvalueB_NB_mean)
      specie2 = specie
      row.names(result)= var
      result2 = cbind(specie2,result)
      result3 = cbind(part,result2)
      result4 = cbind(var,result3)
      log2foldchange_Medium <- rbind(log2foldchange_Medium,result4)
    }}}

nrow(log2foldchange_Medium)
Nutrition = data.frame(rep('Medium', 1400))
log2foldchange_Medium2 = data.frame(log2foldchange_Medium,Nutrition)
colnames(log2foldchange_Medium2) = c('ASV','Compartmert', 'species', 'diff', 'pvalue', 'Nutrition')


log2foldchange_Low <- c()
for(part in melted_tableSyncom_Low[melted_tableSyncom_Low$Species!='LA7339',]$Compartment  %>% unique){
  melted_sub <- melted_tableSyncom_Low[melted_tableSyncom_Low$Species!='LA7339',] %>% subset(Compartment  ==  part) %>% droplevels
  for(specie in melted_sub$Species  %>% unique){
    melted_sub2 <- melted_sub %>% subset(Species  ==  specie) %>% droplevels
    for(var in melted_sub2$variable  %>% unique){
        melted_sub3 <- melted_sub2 %>% subset(variable ==  var) %>% droplevels #'ASV522'
        
        NB = data.frame(t(melted_sub3[melted_sub3$Nutrient  =='High',][,7]))
        B = data.frame(t(melted_sub3[melted_sub3$Nutrient  =='Low',][,7]))
        
        B_mean = apply(B, 1, mean) 
        NB_mean = apply(NB, 1, mean) 
        
        B_NB_mean <- log2(B_mean+1) - log2(NB_mean+1) 
        
        B_NB_mean_statistic = t.test(B, NB)
        
        pvalueB_NB_mean = B_NB_mean_statistic$p.value
        
        result = cbind(B_NB_mean,pvalueB_NB_mean)
        specie2 = specie
        row.names(result)= var
        result2 = cbind(specie2,result)
        result3 = cbind(part,result2)
        result4 = cbind(var,result3)
        log2foldchange_Low <- rbind(log2foldchange_Low,result4)
      }}}

log2foldchange_Low

log2foldchange_LowLA <- c()
for(part in melted_tableSyncom_Low[melted_tableSyncom_Low$Species=='LA7339',]$Compartment  %>% unique){
  melted_sub <- melted_tableSyncom_Low[melted_tableSyncom_Low$Species=='LA7339',] %>% subset(Compartment  ==  part) %>% droplevels
  for(specie in melted_sub$Species  %>% unique){
    melted_sub2 <- melted_sub %>% subset(Species  ==  specie) %>% droplevels
    for(var in melted_sub2$variable  %>% unique){
      melted_sub3 <- melted_sub2 %>% subset(variable ==  var) %>% droplevels #'ASV522'
      
      NB = data.frame(t(melted_sub3[melted_sub3$Nutrient  =='High',][,7]))
      B = data.frame(t(melted_sub3[melted_sub3$Nutrient  =='Low',][,7]))
      
      B_mean = apply(B, 1, mean) 
      NB_mean = apply(NB, 1, mean) 
      
      B_NB_mean <- log2(B_mean+1) - log2(NB_mean+1) 
      
      #B_NB_mean_statistic = t.test(B, NB)
      
      pvalueB_NB_mean = 1
      
      result = cbind(B_NB_mean,pvalueB_NB_mean)
      specie2 = specie
      row.names(result)= var
      result2 = cbind(specie2,result)
      result3 = cbind(part,result2)
      result4 = cbind(var,result3)
      log2foldchange_LowLA <- rbind(log2foldchange_LowLA,result4)
    }}}

log2foldchange_Low2 =rbind(log2foldchange_Low,log2foldchange_LowLA)

nrow(log2foldchange_Low2)
Nutrition = data.frame(rep('Low', 1400))
log2foldchange_Low3 = data.frame(log2foldchange_Low2,Nutrition)
colnames(log2foldchange_Low3) = c('ASV','Compartmert', 'species', 'diff', 'pvalue', 'Nutrition')

log2foldchange = rbind(log2foldchange_Medium2,log2foldchange_Low3)
log2foldchange$species
log2foldchange$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange$Significance[which(log2foldchange$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange$Significance <- log2foldchange$Significance %>% factor
log2foldchange$diff = as.numeric(log2foldchange$diff)

head(log2foldchange)
#Aggregate the dataframe to display as heatmap
display <- log2foldchange  %>%
  acast(formula = species~ASV, mean,
        value.var = "diff")

display2 = data.frame(t(display))


# Obtain the dendrogram
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



log2foldchange$ASV = factor(log2foldchange$ASV, c(ASV_order))
log2foldchange$species = factor(log2foldchange$species, c( "SP9509","SI7820","LP0049", "LJ9250","LP8539", "LA7339","PS"))
log2foldchange$Compartmert = factor(log2foldchange$Compartmert, c('Root','Front'))
log2foldchange$Nutrition = factor(log2foldchange$Nutrition, c('Medium','Low'))

plt_hmap = ggplot(data = log2foldchange[log2foldchange$species!='PS',], aes(species,ASV)) + 
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.9,height = 0.95) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Compartmert*Nutrition, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-11,11),na.value = "#D9D9D9",name = "Fold Change") +
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


