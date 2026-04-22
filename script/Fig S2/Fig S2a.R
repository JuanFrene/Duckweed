library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(pheatmap)
library(reshape2)
library(dplyr)

setwd("C:/Users/juanp/Google Drive/labs/Nottingham/Duckweed/Root traits/Exp1/")

TablaCompleta = read.table("Root traits.txt", header = TRUE, row.names = 1)
colnames(TablaCompleta)

TablaCompleta$Microbiome[which(TablaCompleta$Microbiome == 'B+')] <- "Natural Water"
TablaCompleta$Microbiome[which(TablaCompleta$Microbiome == 'NB')] <- "Sterile Water"
colnames(TablaCompleta) = c("Species","Microbiome","Subsample","CCL","Stele diameter","Endodermis diameter" ,"CCL1 diameter",       
                            "CCL2 diameter","CCL3 diameter","CCL4 diameter","CCLe1 diameter","CCLe2 diameter", "Root diameter",        
                            "Stele N cell","Endodermis N cell","CCL1 N cell","CCL2 N cell","CCL3 N cell","CCL4 N cell",
                            "CCLe1 N cell","CCLe2 N cell","Aerenchyma-like mean area","Aerenchyma-like N","N roots")

TablaCompleta_mean = TablaCompleta%>%group_by(Species,Microbiome)%>%
  summarise_all(mean)

#Table No bacteria sterile water
TablaCompleta_mean_NB = TablaCompleta_mean[TablaCompleta_mean$Microbiome!="Natural Water",]
# Desired order

# Apply factor order
TablaCompleta_mean_NB <- TablaCompleta_mean_NB %>%
  arrange(desc(CCL))

#Order
Species_order = TablaCompleta_mean_NB$Species

# Create matrix after ordering
maxTab2 <- as.matrix(
  log10(TablaCompleta_mean_NB[, c(4:12,14:24)] + 1),
  rownames.force = TRUE
)

# Annotation after ordering
annotation <- data.frame(Species = TablaCompleta_mean_NB$Species)
rownames(annotation) <- colnames(maxTab2)

pheatmap(t(normalize(maxTab2)), 
         cutree_rows = 3, 
         border_color = "grey60",
         cluster_cols = F,
         legend = TRUE,
         size = 0.5,
         annotation_names_col = TRUE,
         annotation_col = annotation,
         #annotation_colors = ann_colors,
         show_colnames = T,
         fontsize_row = 5,fontsize_col = 5)

dev.off()

#Table with bacteria natural water
TablaCompleta_mean_B = TablaCompleta_mean[TablaCompleta_mean$Microbiome=="Natural Water",]
# Desired order

# Apply factor order
TablaCompleta_mean_B$Species <- factor(TablaCompleta_mean_B$Species,
                                       levels = Species_order)

# Reorder dataframe
TablaCompleta_mean_B <- TablaCompleta_mean_B %>%
  arrange(Species)

# Create matrix after ordering
maxTab <- as.matrix(
  log10(TablaCompleta_mean_B[, c(4:12,14:24)] + 1),
  rownames.force = TRUE
)

# Annotation after ordering
annotation <- data.frame(Species = TablaCompleta_mean_B$Species)
rownames(annotation) <- colnames(maxTab2)

pheatmap(t(normalize(maxTab)), 
         cutree_rows = 3, 
         border_color = "grey60",
         cluster_cols = F,
         legend = TRUE,
         size = 0.5,
         annotation_names_col = TRUE,
         annotation_col = annotation,
         #annotation_colors = ann_colors,
         show_colnames = T,
         fontsize_row = 5,fontsize_col = 5)


