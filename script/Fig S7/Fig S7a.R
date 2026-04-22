packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest",'dada2', "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", "reshape2")
sapply(packages, require, character.only = TRUE)              

setwd("G:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/")

# 1. Taxonomy Table
taxa <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/tax_final250_220.rds")
# check your taxonomic classifications #
taxa.print<- taxa
rownames(taxa.print)
head(taxa.print)
unname(head(taxa))

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

###ASV table G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/
asv.table<- readRDS("G:/My Drive/labs/Nottingham/Duckweed/16S Exp1/seqtab_final250_220.rds")
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

#Now we can make the phyloseq object
ps2 <- phyloseq(asv.table2, tax_table(taxa), sample_data(meta))

# while it doesn't seem to be the case, if it were we would use the following lines to remove samples with low sequence counts
set.seed(500)
ps.2 = subset_taxa(ps2, Kingdom == "Bacteria" | Kingdom == "Archaea")
ps.pruned <- prune_samples(sample_sums(ps.2)>=2, ps.2)
ps.3 = subset_samples(ps.pruned, Compartment !=  "BLANK")
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))

ps.perc.3 <- transform_sample_counts(ps.3.SC, function(x) x / sum(x)) 


#PS,LY5280,SI9227,SP7498,SP9509,SI7820,SP5543,SP9192,LM7200,LT9243,LP7760,LV7650,LJ9250,LP0049,LM8389,LT9109,LA7339,LM7123,LP8539
ps.perc.3.LP85392 = subset_samples(ps.perc.3, Specie.1 ==  "LP8539")

#####PCoA for ASV-level data with Bray-Curtis
PCoA <- ordinate(ps.perc.3.LP85392, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x,ps.perc.3.H@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
nrow(x)
'bronw'
cols <- c("Frond" = "#117733", "Root" = "#8B4513", "Water" = "#87CEFA", "Inoculum" = "black")

x$Compartment[which(x$Compartment == 'Front')] <- "Frond"
x2 = x[x$Compartment!='Water',]
unique(x$Specie.1)



LP8539 = ggplot(x, aes(Axis.1, Axis.2)) + #
  geom_point(aes(color=Compartment), size=3) + #
  theme_few() +
  theme(legend.position="top")+ 
  scale_colour_manual(values = cols) +
  labs(y=ylab_text, x=xlab_text)


ggarrange(PS,LY5280,SI9227,SP7498,SP9509,
          SI7820,SP5543,SP9192,LM7200,LT9243,
          LP7760,LV7650,LJ9250,LP0049,LM8389,LT9109,
          LA7339,LM7123,LP8539, 
          ncol = 5, nrow=4, title = 'Growth curves',
          font.label = list(size = 6, face = "bold", color ="red"))






ps.perc.3.LP8539 = subset_samples(ps.perc.3, Specie.1 ==  "LP8539")

# Calculate bray curtis distance matrix
ps.3_bray <- phyloseq::distance(ps.perc.3.LP8539, method = "bray")
# make a data frame from the sample_data
sample_ps.3c <- data.frame(sample_data(ps.perc.3))
# Adonis test
adonis2(ps.3_bray ~ Compartment , data = sample_ps.3c)
#PS: R = 0.159, P = 0.193
#LY5280: R = 0.144, P = 0.633
#SI9227: R = 0.18754, P = 0.385
#SP7498: R = 0.23, P = 0.01
#SP9509: R = 0.159, P = 0.23

#SI7820: R = 0.191, P= 0.027
#SP5543: R = 0.198, P= 0.064
#SP9192: R = 0.293, P= 0.006
#LM7200: R = 0.161, P = 0.342
#LT9243: R = 0.209, P = 0.079

#LP7760: R = 0.198, P = 0.045
#LV7650: R = 0.181, P = 0.114
#LJ9250: R = 0.204, P = 0.007
#LP0049: R = 0.193, P = 0.078
#LM8389: R = 0.072, P = 0.072

#LT9109: R = 0.132, P = 0.463
#LA7339: R = 0.178, P = 0.185
#LM7123: R = 0.197, P = 0.031
#LP8539: R = 0.214, R = 0.076