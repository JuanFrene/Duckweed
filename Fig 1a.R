# install necessary packages 
install.packages("remotes")
remotes::install_github("GuangchuangYu/treeio")
install.packages("BiocManager")
BiocManager::install("ggtree")
saBiocManager::install("DECIPHER")
install.packages("ape")
install.packages('seqinr')
install.packages('adegenet')
install.packages('viridis')
install.packages('ggplot2')
install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.3.tar.gz", repos = NULL)
remotes::install_github("YuLab-SMU/tidytree")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtreeExtra")

install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.3.tar.gz", repos = NULL)
#quit and restart your R and use 
BiocManager::install("ggtree",force=TRUE)

install.packages("dplyr") # To manipulate dataframes
install.packages("readxl") # To read Excel files into R
install.packages("GenomeInfoDb") # for high quality graphics
install.packages("rmdformats")
install.packages("kableExtra")
install.packages("gplots")

#latest development version of "dplyr" from github
if (packageVersion("devtools") < 1.6) {install.packages("devtools")}
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")

# Install phyloseq
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")


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
setwd("G:/My Drive/labs/Nottingham/Duckweed/rpl16 sequence/")
getwd()
register(SnowParam(6))




##################################################################################################################################################
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet("all the sequences.fas", format = "fasta")

# look at some of the sequences (optional)
seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
#BrowseSeqs(aligned, highlight=2)

# write the alignment to a new FASTA file
writeXStringSet(aligned, file="Duckweed_onlyspe_fasta_aligned.fa")

# read in the aligned data
dna <- read.alignment("Duckweed_onlyspe_fasta_aligned.fa", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(dna, matrix = "similarity")

#Write distance table to file and save
df <- melt(as.matrix(D), varnames = c("row", "col"))
df[df$row > df$col,]
#write.csv(df, file ="distance_table.csv")
#grep 'NA' distance_table.csv | cut -d, -f1 --complement | cut -d, -f3 --complement | sed 's/","/ /g' | sed 's/"//g' | tr -s ' ' '\n' | sort | uniq -c -d > problem_fasta_sequences

#View distance matrix. Darker shades of gray mean a larger distance # you can also make cool color plots but they're much more complicated because they use the image() function
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5) + scale_color_viridis()

#create tree object
treeD <- nj(D)
class(treeD) #all trees created using {ape} package will be of class phylo
treeD <- root(treeD, out = 0.1)

treeD <- ladderize(treeD)

#check the formatting of the tip.labels in the tree
treeD$tip.label

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ metadata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import metadata from excel file- Linux
metadata <- read.csv("metadata.csv", sep=',')
# inspect the first 6 Sample_IDs
head(metadata$Sample_ID) 

#determine if all samples are present in the tree file and vice versa
metadata$Sample_ID %in% treeD$tip.label
treeD$tip.label %in% metadata$Sample_ID

#show any sample IDs that are not on the tree
metadata$Sample_ID[!treeD$tip.label %in% metadata$Sample_ID]#character(0) = none


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mypal=c('#d4c29e','#de5275','#6aceee','#8e67a0','#e9ae8c',"#55003b",'#a0ae00','#e0e004','#195850','#535121','#e0aa02','#f35828',"#c1204c",'#5e3a49',"black","white")

#linear with sequence lengths
dat <- data.frame(id=c(1, 25,23,4,5), type=c("Pistia",'Lemna','Spirodela','Landoltia','Landoltia'))

f = ggtree(treeD, layout="circular", color='black', size=0.3) %<+% metadata + 
  geom_tiplab(aes(label=Species), offset=0.01, hjust = 0.5, vjust=0.5, size=3.5, align = TRUE, linetype = "dotted", linesize = 0.15) +
  scale_fill_manual(values = mypal, guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))+ 
  geom_treescale(fontsize=3, linesize=0.5)+
  geom_hilight(data=dat, mapping=aes(node=id, fill=type),
               align="right",extend=0.3,show.legend=TRUE,alpha=0.3) 
dat3
dat3 = metadata
colnames(dat3) = c('Sample_ID','Species','Species_ID','Root_complexity2','Nroot2')

g = f + geom_fruit(data=dat3, geom=geom_bar,
               mapping=aes(y=Sample_ID, x=Nroot2),
               pwidth=0.5, 
               offset = 0.7,
               orientation="y", 
               stat="identity",
               axis.params=list(
                 axis       = "x",
                 text.size  = 1.8,
                 hjust      = 1,
                 vjust      = 0.5,
                 nbreak     = 8),
               grid.params=list()) 
g
dat4 = metadata
colnames(dat4) = c('Sample_ID','Species','Species_ID','Root_complexity2','Nroot2')

g + geom_fruit(data=dat4, geom=geom_bar,
                   mapping=aes(y=Sample_ID, x=Root_complexity2),
                   pwidth=0.5, 
                   offset = 1,
                   orientation="y", 
                   stat="identity",
                   axis.params=list(
                     axis       = "x",
                     text.size  = 1,5,
                     hjust      = 1,
                     vjust      = 0.5,
                     nbreak     = 6),
                   grid.params=list()) 



f = ggtree(treeD, color='black', size=0.2, ladderize = FALSE) %<+% metadata + 
  geom_tiplab(aes(label=Species_ID), offset=0.01, hjust = 0.5, vjust=0.5, size=3.5, align = FALSE, linetype = "dotted", linesize = 0.15) +
  scale_fill_manual(values = mypal, guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))+ 
  theme_tree2()+
  geom_highlight(data=dat, 
                 mapping=aes(node=id, fill=type),
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) 
