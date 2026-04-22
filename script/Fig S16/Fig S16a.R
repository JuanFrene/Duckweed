packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'dada2')
sapply(packages, require, character.only = TRUE)              

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/DExp2 16S analysis/") #Users/juanp/Google

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

ps.Syncom.SC = subset_samples(ps.pruned, Species !=  "Control")
ps.Syncom.SC2 = subset_samples(ps.Syncom.SC, Species !=  "Inoculum")

ps.Syncom.Syncom = subset_samples(ps.Syncom.SC2, Treatment ==  "SynCom")
ps.Syncom.Syncom2 <- prune_taxa(taxa_sums(ps.Syncom.Syncom)>1, ps.Syncom.Syncom)

ps.Syncom.Syncom2 <- prune_taxa(taxa_sums(ps.Syncom.Syncom)>1, ps.Syncom.Syncom)
ps.Syncom.Syncom3 <- prune_samples(sample_sums(ps.Syncom.Syncom2)>1, ps.Syncom.Syncom2) 

ps.Syncom.Front = subset_samples(ps.Syncom.Syncom, Compartment ==  "Front")
ps.Syncom.Front2 <- prune_taxa(taxa_sums(ps.Syncom.Front)>0, ps.Syncom.Front)
ps.Syncom.Front3 <- prune_samples(sample_sums(ps.Syncom.Front)>0, ps.Syncom.Front) 
ps.Frond.perc <- transform_sample_counts(ps.Syncom.Front3, function(x) x / sum(x)) 


ps.Syncom.Root = subset_samples(ps.Syncom.Syncom, Compartment ==  "Root")
ps.Syncom.Root2 <- prune_taxa(taxa_sums(ps.Syncom.Root)>0, ps.Syncom.Root)
ps.Syncom.Root3 <- prune_samples(sample_sums(ps.Syncom.Root2)>0, ps.Syncom.Root2) 

ps.Syncom.Root4 <- subset_samples(ps.Syncom.Root3, Species != 'PS')
ps.Root.perc <- transform_sample_counts(ps.Syncom.Root4, function(x) x / sum(x)) 

###CAP
data = ps.Syncom.Root3@sam_data

CAP2 <- ordinate(ps.Root.perc, "CAP", "bray",~data$Nutrient)
summary(CAP2)
xxlab_text = 'CAP1 (46.63%)'
ylab_text = 'CAP2 (23.86%)'

x<-data.frame(CAP2$CCA$wa[, 1:2])
x <- merge(x,ps.Syncom.Root3@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
nrow(x)

safe_colorblind_palette <- c("#44AA99", "#999933", "#882255", "#332288")

safe_colorblind_palette2 <- c("#1912f1", "#117733","#DDCC77","#882615" )

x$Complexity <- "High"
x$Complexity[which(x$Species == 'LA7339')] <- "Low"
x$Complexity[which(x$Species == 'LP8539')] <- "Low"
x$Complexity[which(x$Species == 'LJ9250')] <- "Medium"
x$Complexity[which(x$Species == 'LP0049')] <- "Medium"
x$Complexity <- x$Complexity %>% factor

ggplot(x, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Nutrient), size=4,pch=21, color="black") +
  theme_few() +
  geom_tile()+
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  #labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette)+
  geom_text(aes(label = "p = 0.776",
    x = -0.6, y = 0.63, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.032",
                x = -0.583, y = 0.57, size = 6, fill = "black"))

ggplot(x, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Complexity ), size=4,pch=21, color="black") +
  theme_few() +
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  #labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette2)+
  geom_text(aes(label = "p = 0.049",
                x = -0.6, y = 0.63, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.038",
                x = -0.587, y = 0.57, size = 6, fill = "black"))



###CAP Frond
dataF = ps.Frond.perc@sam_data
CAPF <- ordinate(ps.Frond.perc, "CAP", "bray",~dataF$Nutrient)
summary(CAPF)
xlab_text = 'CAP1 (33,06%)'
ylab_text = 'CAP2 (11.49%)'

f<-data.frame(CAPF$CCA$wa[, 1:2])
f <- merge(f,ps.Syncom.Front3@sam_data, by=0)
rownames(f)<-x$Row.names
f$Row.names <- NULL
nrow(f)

f$Complexity <- "High"
f$Complexity[which(f$Species == 'LA7339')] <- "Low"
f$Complexity[which(f$Species == 'LP8539')] <- "Low"
f$Complexity[which(f$Species == 'LJ9250')] <- "Medium"
f$Complexity[which(f$Species == 'LP0049')] <- "Medium"
f$Complexity <- f$Complexity %>% factor

safe_colorblind_palette <- c("#44AA99", "#999933", "#882615", "#332288", "#6699CC", "#888888")

ggplot(f, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Nutrient), size=4,pch=21, color="black") +
  theme_few() +
  geom_tile()+
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette)+
  geom_text(aes(label = "p = 0.008",
                x = -0.443, y = 0.66, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.044",
                x = -0.43, y = 0.6, size = 6, fill = "black"))

ggplot(f, aes(CAP1, CAP2)) + #
  geom_point(aes(fill=Complexity ), size=4,pch=21, color="black") +
  theme_few() +
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  labs(y=ylab_text, x=xlab_text)+
  scale_fill_manual(values = safe_colorblind_palette2)+
  geom_text(aes(label = "p = 0.297",
                x = -0.443, y = 0.66, size = 6, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R2 = 0.013",
                x = -0.433, y = 0.6, size = 6, fill = "black"))



