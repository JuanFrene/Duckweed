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

####Total_vel 
log2foldchange_TotalVel <- c()
for(specie in Model$Species  %>% unique){
  melted_sub <- Model %>% subset(Species ==  specie) %>% droplevels
  
  NB = data.frame(t(melted_sub[melted_sub$Treatment=='No-Microbiome',][,19]))
  B = data.frame(t(melted_sub[melted_sub$Treatment=='With-Microbiome',][,19]))
  
  NB_mean = apply(NB, 1, mean) 
  B_mean = apply(B, 1, mean) 
  
  Fmax <- log2(B_mean) - log2(NB_mean) 
  
  Fmax_statistic = t.test(B,NB)
  Fmax_pvalue = Fmax_statistic$p.value#Water_mean+1
  
  result = cbind(Fmax,Fmax_pvalue)
  row.names(result)= specie
  log2foldchange_TotalVel <- rbind(log2foldchange_TotalVel,result)
}

data.frame(log2foldchange_TotalVel)

log2foldchange_TotalVel = cbind(row.names(log2foldchange_TotalVel),log2foldchange_TotalVel)
colnames(log2foldchange_TotalVel) = c('Species','Fmax','Fmax_pvalue')

Model_TotalVel = merge(Model[,c(2,3,19)],log2foldchange_TotalVel, by='Species')

pval_thres <- 0.05
Model_TotalVel$Significance <- "No Significant"
Model_TotalVel$Significance[which(Model_TotalVel$Fmax_pvalue < pval_thres)] <- "Significant"
Model_TotalVel$Significance <- Model_TotalVel$Significance %>% factor
head(Model_TotalVel)

Model_TotalVel$Log <- log10(Model_TotalVel$Total_vel+1)

Model_TotalVel$Treatment[Model_TotalVel$Treatment=="No-Microbiome"]<-"Sterile Water"
Model_TotalVel$Treatment[Model_TotalVel$Treatment=="With-Microbiome"]<-"Natural Water"

Model_TotalVel$Treatment <- factor(Model_TotalVel$Treatment,c('Sterile Water','Natural Water'))


ggplot(data = Model_TotalVel, aes(Treatment,Species)) +
  geom_raster(aes(fill = Log)) +
  geom_tile(aes(color = Treatment),fill = 'transparent', linewidth = 0.2,width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  scale_color_manual(values = c('black',"black"), aesthetics = "colour",name = "Significance vs Sterile Water        ") + #Significance Genotype vs Col-0
  theme(size_border = 2.5) +
  theme_few()

ncol(Model)
####Initial_vel
log2foldchange_InitialVel <- c()
for(specie in Model$Species  %>% unique){
  melted_sub <- Model %>% subset(Species ==  specie) %>% droplevels
  
  NB = data.frame(t(melted_sub[melted_sub$Treatment=='No-Microbiome',][,20]))
  B = data.frame(t(melted_sub[melted_sub$Treatment=='With-Microbiome',][,20]))
  
  NB_mean = apply(NB, 1, mean) 
  B_mean = apply(B, 1, mean) 
  
  InitialVel <- log2(B_mean) - log2(NB_mean) 
  
  InitialVel_statistic = t.test(B,NB)
  InitialVel_pvalue = InitialVel_statistic$p.value#Water_mean+1
  
  result = cbind(InitialVel,InitialVel_pvalue)
  row.names(result)= specie
  log2foldchange_InitialVel <- rbind(log2foldchange_InitialVel,result)
}

data.frame(log2foldchange_InitialVel)

log2foldchange_InitialVel = cbind(row.names(log2foldchange_InitialVel),log2foldchange_InitialVel)
colnames(log2foldchange_InitialVel) = c('Species','InitialVel','InitialVel_pvalue')

Model_InitialVel = merge(Model[,c(2,3,20)],log2foldchange_InitialVel, by='Species')

pval_thres <- 0.05
Model_InitialVel$Significance <- "No Significant"
Model_InitialVel$Significance[which(Model_InitialVel$InitialVel_pvalue < pval_thres)] <- "Significant"
Model_InitialVel$Significance <- Model_InitialVel$Significance %>% factor
head(Model_InitialVel)

Model_InitialVel$Log <- log10(Model_InitialVel$Initial_vel+1)

Model_InitialVel$Treatment[Model_InitialVel$Treatment=="No-Microbiome"]<-"Sterile Water"
Model_InitialVel$Treatment[Model_InitialVel$Treatment=="With-Microbiome"]<-"Natural Water"

Model_InitialVel$Treatment <- factor(Model_InitialVel$Treatment,c('Sterile Water','Natural Water'))

ggplot(data = Model_InitialVel, aes(Treatment,Species)) +
  geom_raster(aes(fill = Log)) +
  geom_tile(aes(color = Treatment),fill = 'transparent', linewidth = 0.2,width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  scale_color_manual(values = c('black',"black"), aesthetics = "colour",name = "Significance vs Sterile Water        ") + #Significance Genotype vs Col-0
  theme(size_border = 2.5) +
  theme_few()


####Tfmax
log2foldchange_Tfmax <- c()
for(specie in Model$Species  %>% unique){
  melted_sub <- Model %>% subset(Species ==  specie) %>% droplevels
  
  NB = data.frame(t(melted_sub[melted_sub$Treatment=='No-Microbiome',][,15]))
  B = data.frame(t(melted_sub[melted_sub$Treatment=='With-Microbiome',][,15]))
  
  NB_mean = apply(NB, 1, mean) 
  B_mean = apply(B, 1, mean) 
  
  Tfmax <- log2(B_mean) - log2(NB_mean) 
  
  Tfmax_statistic = t.test(B,NB)
  Tfmax_pvalue = Tfmax_statistic$p.value#Water_mean+1
  
  result = cbind(Tfmax,Tfmax_pvalue)
  row.names(result)= specie
  log2foldchange_Tfmax <- rbind(log2foldchange_Tfmax,result)
}

data.frame(log2foldchange_Tfmax)

log2foldchange_Tfmax = cbind(row.names(log2foldchange_Tfmax),log2foldchange_Tfmax)
colnames(log2foldchange_Tfmax) = c('Species','Tfmax','Tfmax_pvalue')

Model_Tfmax = merge(Model[,c(2,3,16)],log2foldchange_Tfmax, by='Species')

pval_thres <- 0.05
Model_Tfmax$Significance <- "No Significant"
Model_Tfmax$Significance[which(Model_Tfmax$Tfmax_pvalue < pval_thres)] <- "Significant"
Model_Tfmax$Significance <- Model_Tfmax$Significance %>% factor
head(Model_Tfmax)

Model_Tfmax$Log <- log10(Model_Tfmax$Tfmax.x)

Model_Tfmax$Treatment[Model_Tfmax$Treatment=="No-Microbiome"]<-"Sterile Water"
Model_Tfmax$Treatment[Model_Tfmax$Treatment=="With-Microbiome"]<-"Natural Water"

Model_Tfmax$Treatment <- factor(Model_Tfmax$Treatment,c('Sterile Water','Natural Water'))


ggplot(data = Model_Tfmax, aes(Treatment,Species)) +
  geom_raster(aes(fill = Tfmax.x)) +
  geom_tile(aes(color = Treatment),fill = 'transparent', linewidth = 0.2,width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_isoluminant_cjm_75_c24") +
  scale_color_manual(values = c('black',"black"), aesthetics = "colour",name = "Significance vs Sterile Water        ") + #Significance Genotype vs Col-0
  theme(size_border = 2.5) +
  theme_few()



#############
correlation
Model ##Growth parameter table
TablaCompleta ##Root trait parameter table

nrow(TablaCompleta)
nrow(Model)
All = cbind(TablaCompleta, Model[,-(2:4)])
ncol(TablaCompleta)

##Corr y p-value
CorrNB<-cor(All[1:95,-c(1:3,25:35)], method = 'pearson')
p.matnB <-cor_pmat(All[1:95,-c(1:3,25:35)], method = 'pearson')

CorrB<-cor(All[-(1:95),-c(1:3,25:35)], method = 'pearson')
p.matB <-cor_pmat(All[-(1:95),-c(1:3,25:35)], method = 'pearson')


in_vel = data.frame(cbind(rownames(CorrB), CorrB[,24], p.matB[,24]))
colnames(in_vel) = c('Root','r', 'p.value')
in_vel$r = as.numeric(in_vel$r)

####Adjust the p-value
in_vel$SignPvalue <- "NoSignificant"
pval_thres <- 0.05
in_vel$SignPvalue[which(in_vel$p.value < pval_thres)] <- "p < 0.05"
in_vel$SignPvalue <- in_vel$SignPvalue %>% factor

ggplot(in_vel[-(22:27),], aes(r,Root)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_few()+ 
  theme(legend.position="top")+
  geom_tile(aes(color = SignPvalue),fill = '#00000000', linewidth = 0.2,width = 0.01,height = 0.95) +
  xlab('r')+
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significative")
  