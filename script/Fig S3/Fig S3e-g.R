library(MicEco)
library(phyloseq)
library(microbiome)

setwd("C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/")

# Taxonomy from dada
taxa <- readRDS("C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/tax_final250_220.rds")
# check your taxonomic classifications #
taxa.print<- taxa
rownames(taxa.print)
head(taxa.print)

###Metadata
meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)

###ASV table from dada
asv.table<- readRDS("C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/16S Exp1/seqtab_final250_220.rds")
asv.table2<- otu_table(asv.table, taxa_are_rows=FALSE)

# 13. Make a Phyloseq object ####
ps.1 <- phyloseq(asv.table2, tax_table(taxa), sample_data(meta))


#Phyloseq selection
ps.2 = subset_taxa(ps.1, Kingdom == "Bacteria" | Kingdom == "Archaea")
ps.pruned <- prune_samples(sample_sums(ps.2)>=2, ps.2)
ps.perc <- transform_sample_counts(ps.pruned, function(x) x / sum(x)) 
ordered(sample_sums(ps.pruned))

#Delete the blank
ps.3 = subset_samples(ps.pruned, Compartment !=  "BLANK")
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))

###Core microbiome by species.
##ps.3.9243
ps.3.perc <- transform_sample_counts(ps.3, function(x) x / sum(x)) 
ps.3.core <- core(ps.3.perc, detection = 1/1000, prevalence = 50/100)

#Diagram venn based on the core
ps_venn(ps.3.core, group = "Compartment", fill = c("red", "blue"), quantities = list(type=c("counts")))
