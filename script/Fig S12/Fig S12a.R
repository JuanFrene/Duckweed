###################################################################################################################################################
#https://yulab-smu.top/treedata-book/
#https://epirhandbook.com/en/phylogenetic-trees-1.html
#https://yulab-smu.top/treedata-book/chapter4.html
#https://yulab-smu.top/treedata-book/chapter10.html
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/advanced-models-for-differential-abundance.html#deseq2-differential-abundance-testing-for-sequencing-data
#https://rstudio.com/products/rstudio/download/#download
#https://cran.r-project.org/bin/windows/Rtools/
#####################################################################################################################################################
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
library("BiocParallel")
library("readxl")
library(ggtree)
library(ggtreeExtra)
library(dplyr)
library(reshape2)
library(treeio)
library(readr)
library(ggstar)
library(Biostrings)
library(ggnewscale)
#options(getClass.msg=FALSE)

###################################################################################################################################################
# Set working directories- Windows-HP
setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/Syncom/")
getwd()
register(SnowParam(6))
##################################################################################################################################################
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet("Syncom_Sequence.txt", format = "fasta")

# look at some of the sequences (optional)
seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=2)

# write the alignment to a new FASTA file
writeXStringSet(aligned, file="Duckweed_Syncom_fasta_aligned.fa")

# read in the aligned data
dna <- read.alignment("Duckweed_Syncom_fasta_aligned.fa", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(dna, matrix = "similarity")

#Write distance table to file and save
df <- melt(as.matrix(D), varnames = c("row", "col"))
df[df$row > df$col,]
write.csv(df, file ="distance_table.csv")
#grep 'NA' distance_table.csv | cut -d, -f1 --complement | cut -d, -f3 --complement | sed 's/","/ /g' | sed 's/"//g' | tr -s ' ' '\n' | sort | uniq -c -d > problem_fasta_sequences

#View distance matrix. Darker shades of gray mean a larger distance # you can also make cool color plots but they're much more complicated because they use the image() function
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5) + scale_color_viridis()

#create tree object
treeD <- nj(D)
class(treeD) #all trees created using {ape} package will be of class phylo

treeD <- ladderize(treeD)

#check the formatting of the tip.labels in the tree
treeD$tip.label

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ metadata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import metadata from excel file- Linux
metadata2 <- read.csv("Duckweed syncom desing.csv", sep=',')
# inspect the first 6 Sample_IDs
head(metadata2$Samples_ID) 
nrow(metadata2)
#determine if all samples are present in the tree file and vice versa
metadata2$Samples_ID %in% treeD$tip.label
treeD$tip.label %in% metadata2$Samples_ID

#show any sample IDs that are not on the tree
metadata2$Samples_ID[!tree$tip.label %in% metadata2$Samples_ID]#character(0) = none
Class_order
unique(metadata2$Class)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mypal=c('#d4c29e','#de5275','#6aceee','#8e67a0','#e9ae8c',"#55003b",'#a0ae00','#e0e004','#195850','#535121','#e0aa02','#f35828',"#c1204c",'#5e3a49',"black","white")


ggtree(treeD, layout="circular", color='black', size=0.25) %<+% metadata2 + 
  geom_tiplab(aes(label=Syncom_ID), offset=0.01, hjust = 0.5, vjust=0.5, size=3.5, align = T, linetype = "dotted", linesize = 0.05) +
  scale_fill_manual(values = mypal, guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))+ 
  geom_treescale(fontsize=3, linesize=0.4)+
  geom_tippoint(aes(colour=Class))
  
