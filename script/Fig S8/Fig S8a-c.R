library(devtools)
library(MicEco)
library(phyloseq)
library(microbiome)
library(eulerr)
library(microbiomeutilities)


# Phyloseq object
seqtab <- readRDS("G:/My Drive/labs/Nottingham/Duckweed/Figuras paper/Clean data/Fig S8/Fig S8 a-c/seqtab Fig S8.rds")

###Core microbiome by species.
table(meta(ps.3.SCSI)$Specie.1, useNA = "always")

ps.3.n = subset_samples(ps.3.SCSI, Specie.1 ==  "SP7498")

uniques_ASV <- unique(as.character(meta(ps.3.n)$Compartment))
print(uniques)

list_core <- c() # an empty object to store information

for (n in uniques_ASV){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.3.n, Compartment == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.00001, # 0.001 in atleast 90% samples 
                         prevalence = 0.01)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

Uniques <- c(Root="#d6e2e9", water="#cbf3f0", Frond="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)
