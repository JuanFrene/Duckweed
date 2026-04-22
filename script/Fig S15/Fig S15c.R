library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(ggdendro)
library(factoextra)
library(dplyr)
library(ggdendro)
library(reshape2)
library(cowplot)
library(pals)
library(paletteer)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/")

Status <- read.table("Exp2_Duckweed_fitness.txt", header = TRUE, row.names = 1)

Status$Species <- factor(Status$Species,rev(c('SP9509','SI7820','LJ9250','LP0049','LA7339','LP8539')))
Status$Fitness <- factor(Status$Fitness,c('Dead','Stressed','Healthy'))
Status$Treatment[Status$Treatment=="NB"]<-"Sterile Water"

Status$Treatment <- factor(Status$Treatment,c('Sterile Water','Syncom'))
Status[Status$Species=='LA7339',]

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

log2foldchange_Table <- c()

for(treat in Status$Treatment  %>% unique){
  melted_sub <- Status %>% subset(Treatment  ==  treat) %>% droplevels
  for(specie in melted_sub$Species  %>% unique){
    melted_sub2 <- melted_sub %>% subset(Species  ==  specie) %>% droplevels
    for(var in melted_sub2$Fitness  %>% unique){
      melted_sub3 <- melted_sub2 %>% subset(Fitness ==  var) %>% droplevels
      for(nutrient in melted_sub3[melted_sub3$Medium !='Dir',]$Medium  %>% unique){
        melted_sub4 <- melted_sub3[melted_sub3$Medium !='Dir',] %>% subset(Medium  ==  nutrient) %>% droplevels
        
        
        NB = data.frame(t(melted_sub3[melted_sub3$Medium =='Dir',][,9]))
        B = data.frame(t(melted_sub4[,9]))
        
        B_mean = apply(B, 1, mean) 
        NB_mean = apply(NB, 1, mean) 
        
        B_NB_mean <- log2(B_mean+1) - log2(NB_mean+1)+0.1 
        
        B_NB_mean_statistic = t.test(B,NB)
        
        pvalueB_NB_mean = B_NB_mean_statistic$p.value
        
        result = cbind(B_NB_mean,pvalueB_NB_mean)
        specie2 = specie
        row.names(result)= var
        result2 = cbind(specie2,result)
        result3 = cbind(nutrient,result2)
        result4 = cbind(treat,result3)
        log2foldchange_Table <- rbind(log2foldchange_Table,result4)
      }}}}

nrow(log2foldchange_Table)

nrow(log2foldchange_Table)

log2foldchange_Table_Complete = data.frame(cbind(row.names(log2foldchange_Table),log2foldchange_Table))
colnames(log2foldchange_Table_Complete)=c('Fitness','Treatment','Nutrition','Species','diff','pvalue')
log2foldchange_Table_Complete4 =log2foldchange_Table_Complete
head(log2foldchange_Table_Complete4)

log2foldchange_Table_Complete4$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange_Table_Complete4$Significance[which(log2foldchange_Table_Complete4$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange_Table_Complete4$Significance <- log2foldchange_Table_Complete4$Significance %>% factor

#log2foldchange_Table_Complete4$var = paste(log2foldchange_Table_Complete4$Root_traits, log2foldchange_Table_Complete4$Nutrition, sep='_') 
log2foldchange_Table_Complete4$diff = as.numeric(log2foldchange_Table_Complete4$diff)

unique(log2foldchange_Table_Complete4$Root_traits)

log2foldchange_Table_Complete4$Species = factor(log2foldchange_Table_Complete4$Species, rev(c('SP9509', 'SI7820', 'LP0049', 'LJ9250', 'LP8539', 'LA7339')))

log2foldchange_Table_Complete4$Fitness = factor(log2foldchange_Table_Complete4$Fitness, c('Healthy', 'Stressed', 'Dead'))
log2foldchange_Table_Complete4$diff = as.numeric(log2foldchange_Table_Complete4$diff)


ggplot(data = log2foldchange_Table_Complete4, aes(Fitness, Species)) +
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.3,width = 0.9,height = 0.95) + #
  facet_grid(~Nutrition*Treatment, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-10,10),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('gray',"black"),na.value =  "transparent",name = "Significance vs 1/2 SH") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05, size = 10, color='black'),
        axis.title = element_blank(), legend.position = 'right')
