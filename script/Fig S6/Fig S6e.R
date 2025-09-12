library(ggplot2)
library(paletteer)
library(reshape2)
library(dplyr)
library(ggthemes)


setwd("G:/My Drive/labs/Nottingham/Duckweed/Ex2/")

Model <- read.table( "Exp2_Duckweed_Curve_Parameters.txt", header = TRUE, row.names = 1)
head(Model)
Model$Species <- factor(Model$Species,rev(c('PS','SP9509','SI7820','LP0049','LJ9250','LA7339','LP8539')))
Model = Model[Model$Treatment != 'NatCom',]
Model = Model[Model$Species != 'PS',]

####Fmax
log2foldchange_Fmax <- c()
for(specie in Model$Species  %>% unique){
  melted_sub <- Model %>% subset(Species ==  specie) %>% droplevels
  
  NB = data.frame(t(melted_sub[melted_sub$Treatment=='NB',][,5]))
  B = data.frame(t(melted_sub[melted_sub$Treatment=='SynCom',][,5]))
  
  NB_mean = apply(NB, 1, mean) 
  B_mean = apply(B, 1, mean) 
  
  Fmax <- log2(B_mean) - log2(NB_mean) 
  
  Fmax_statistic = t.test(B,NB)
  Fmax_pvalue = Fmax_statistic$p.value#Water_mean+1
  
  result = cbind(Fmax,Fmax_pvalue)
  row.names(result)= specie
  log2foldchange_Fmax <- rbind(log2foldchange_Fmax,result)
}

data.frame(log2foldchange_Fmax)

log2foldchange_Fmax = cbind(row.names(log2foldchange_Fmax),log2foldchange_Fmax)
colnames(log2foldchange_Fmax) = c('Species','Fmax','Fmax_pvalue')

Model_Fmax = merge(Model[,c(1,2,3,5)],log2foldchange_Fmax, by='Species')
colnames(Model_Fmax) = c('Species','Nutient', 'Treatment','Fmax','foldchange','Fmax_pvalue')

pval_thres <- 0.05
Model_Fmax$Significance <- "No Significant"
Model_Fmax$Significance[which(Model_Fmax$Fmax_pvalue < pval_thres)] <- "Significant"
Model_Fmax$Significance <- Model_Fmax$Significance %>% factor
head(Model_Fmax)

Model_Fmax$Log <- log10(Model_Fmax$Max)

Model_Fmax$Treatment[Model_Fmax$Treatment=="NB"]<-"Sterile Water"

Model_Fmax$Treatment <- factor(Model_Fmax$Treatment,c('Sterile Water','SynCom'))


ggplot(data = Model_Fmax, aes(Treatment,Species)) +
  geom_raster(aes(fill = Fmax)) +
  facet_grid(.~ Nutient)+
  geom_tile(aes(color = Treatment),fill = 'transparent', linewidth = 0.2,width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(0,5),na.value = "#D9D9D9",name = "Fmax") +
  scale_color_manual(values = c('black',"black"), aesthetics = "colour",name = "Significance vs Sterile Water        ") + #Significance Genotype vs Col-0
  theme(size_border = 2.5) +
  theme_few()+
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05))

