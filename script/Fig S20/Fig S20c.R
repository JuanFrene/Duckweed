######Load packages
library(Rmisc)
library(reshape2)
library(ggrepel)
library(scales)
library(ggtree)
library(harrietr)
library(emmeans)
library(paletteer)
library(pals)
library(heatmaply)
library(ggcorrplot)
library("viridis")
library(dendextend)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(ggdendro)



######Set work directory
setwd("G:/My Drive/labs/Nottingham/Duckweed/Drop out/Analysis files/")

######Open data file
Data <- read.table("LA7339 dropout.txt", header = TRUE, row.names = 1)
head(Data)
ID=rownames(Data)
Data= cbind(ID,Data)
Data$Syncom[which(Data$Syncom == 'Full_Syncom')] <- "Syncom"


### Now we need to summarize 
melted <- Data[,-(2:11)] %>% melt #Tab_f 
colnames(melted) <- c("ID","parameter","value")

###melted data and metadata
melted <- merge(Data[,c(1,4,10:11)],melted,  by='ID')

#Do a model per ion
Res_em <- NULL
Res_p <- NULL

for(batch in melted$Batch %>% unique){
  melted_batch <- melted %>% subset(Batch  == batch) %>% droplevels
  
  for(ion in melted_batch$variable %>% unique){
    melted_sub <- melted_batch %>% subset(variable  == ion) %>% droplevels
    m4 <- lm(data = melted_sub,
             formula = value ~ Syncom)
    m4_res <- emmeans(m4,pairwise~Syncom,ref =4, adjust = "none")
    m4_em <- m4_res$emmeans %>% as.data.frame
    m4_p <- m4_res$contrasts %>% as.data.frame
    m4_em$Ion <- ion
    m4_p$Ion <- ion
    Res_em <- rbind(Res_em,m4_em)
    Res_p <- rbind(Res_p,m4_p)
  }
}


#Arrange the pvalues dataframe
Res_p$contrast <- Res_p$contrast %>% gsub(pattern = " ",replacement = "")

#Now the contrast for Drop temperature
c9 <- Res_p$contrast %>% grep(pattern = "Syncom",value = T) %>% unique


####Subset the Drop temperature contrast
Res_p_Drop <- which(Res_p$contrast %in% c9) %>%
  Res_p[.,] %>% droplevels


####Adjust the p-value
Res_p_Drop$Significance <- "NoSignificant"
pval_thres <- 0.1
Res_p_Drop$Significance[which(Res_p_Drop$p.value < pval_thres)] <- "q < 0.1"
Res_p_Drop$Significance <- Res_p_Drop$Significance %>% factor

df_temp_Drop <- Res_p_Drop$contrast %>% gsub(pattern = "Complete-",replacement = "") %>%
  gsub(pattern = "()", replacement = "") %>%
  as.data.frame

df_temp_Drop

df_temp_Drop$Syncom <- df_temp_Drop$Syncom %>% gsub(pattern = "-Syncom",replacement = "")


colnames(df_temp_Drop) <- c("Syncom")

Res_p_Drop <- cbind(Res_p_Drop,df_temp_Drop)

Res_p_Drop$UId <- paste(Res_p_Drop$Genotype,Res_p_Drop$Ion,sep = "_")


df_temp_Drop <- data.frame(Syncom = Res_p_Drop$Syncom,
                           p.value =Res_p_Drop$p.value, 
                           estimate = Res_p_Drop$estimate)


df_temp_Drop$Significance <- "NoSignificant"
df_temp_Drop$Significance[which(df_temp_Drop$p.value < pval_thres)] <- "q < 0.1"
df_temp_Drop$Significance <- df_temp_Drop$Significance %>% factor

colnames(df_temp_Drop)='Syncom'

unique(Res_p_Drop$Ion)

Res_p_Drop$Syncom = factor(Res_p_Drop$Syncom, c('R289-Syncom','R245-Syncom',"NB-Syncom","P9-Syncom",'P153-Syncom','P234-Syncom',"P129-Syncom","P118A-Syncom","P162-Syncom","R34-Syncom","R186-Syncom","P128A.22-Syncom","P74A-Syncom","R148A-Syncom","R345.2-Syncom","R151-Syncom","R28-Syncom","P69-Syncom","P192-Syncom","R88-Syncom","P65-Syncom"))

ggplot(data = Res_p_Drop[,-10], aes(Syncom,Ion)) +
  geom_raster(aes(fill = estimate)) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.5,width = 0.85,height = 0.85) + 
  #geom_text(aes(label = Label),color = "black",size = 5)+
  #facet_grid(~Batch,space = "free",scales = "free") +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2.5,2.5),na.value = "#D9D9D9",oob = squish,name = "Abundance \nz-score") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs. Complete Syncom: ") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.text.x = element_blank()#,
        #legend.text = element_text(size = 5)
  ) +
  scale_y_discrete(expand = c(0,0))
