library(ggplot2)
library(paletteer)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp1/")

Model <- read.table( "Growth and vel.txt", header = TRUE, row.names = 1)
head(Model)
Model$Species <- factor(Model$Species,rev(c('PS','SP7820','SI9227','SP9509','SP5543',
                                            'SP7498','SP9192','LP7760','LP0049','LT9109',
                                            'LJ9250','LM7123','LM8389','LV7650','LY5280',
                                            'LM7200','LA7339','LP8539','LT9243')))

####Fmax
log2foldchange_Fmax <- c()
for(specie in Model$Species  %>% unique){
  melted_sub <- Model %>% subset(Species ==  specie) %>% droplevels
  
  NB = data.frame(t(melted_sub[melted_sub$Treatment=='No-Microbiome',][,15]))
  B = data.frame(t(melted_sub[melted_sub$Treatment=='With-Microbiome',][,15]))
  
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

Model_Fmax = merge(Model[,c(2,3,15)],log2foldchange_Fmax, by='Species')

pval_thres <- 0.05
Model_Fmax$Significance <- "No Significant"
Model_Fmax$Significance[which(Model_Fmax$Fmax_pvalue < pval_thres)] <- "Significant"
Model_Fmax$Significance <- Model_Fmax$Significance %>% factor
head(Model_Fmax)

Model_Fmax$Log <- log10(Model_Fmax$Max)

Model_Fmax$Treatment[Model_Fmax$Treatment=="No-Microbiome"]<-"Sterile Water"
Model_Fmax$Treatment[Model_Fmax$Treatment=="With-Microbiome"]<-"Natural Water"

Model_Fmax$Treatment <- factor(Model_Fmax$Treatment,c('Sterile Water','Natural Water'))


ggplot(data = Model_Fmax, aes(Treatment,Species)) +
  geom_raster(aes(fill = Log)) +
  geom_tile(aes(color = Treatment),fill = 'transparent', linewidth = 0.2,width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  scale_color_manual(values = c('black',"black"), aesthetics = "colour",name = "Significance vs Sterile Water        ") + #Significance Genotype vs Col-0
  theme(size_border = 2.5) +
  theme_few()
