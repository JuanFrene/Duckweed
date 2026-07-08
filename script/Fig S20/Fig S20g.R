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


######Set work directory
setwd("G:/My Drive/labs/Nottingham/Duckweed/Drop out/Analysis files/")

######Open data file
Data2 <- read.table("LA7339 mono.txt", header = TRUE, row.names = 1)
head(Data2)
ID=rownames(Data2)
Data2= cbind(ID,Data2)
Data2$Syncom[which(Data2$Syncom == 'Full_Syncom')] <- "Syncom"

### Now we need to summarize 
melted <- Data2 %>% melt  
colnames(melted) <- c("ID",'Syncom',"parameter","value")

###melted data and metadata
#melted <- merge(Data[,1:10],melted,  by='ID')

#Do a model per Root_trait
Res_em <- NULL
Res_p <- NULL

for(Root_trait in melted$parameter %>% unique){
  melted_sub <- melted %>% subset(parameter  == Root_trait) %>% droplevels
  m4 <- lm(data = melted_sub,
           formula = value ~ Syncom)
  m4_res <- emmeans(m4,pairwise~Syncom,ref =4, adjust = "none")
  m4_em <- m4_res$emmeans %>% as.data.frame
  m4_p <- m4_res$contrasts %>% as.data.frame
  m4_em$Root_trait <- Root_trait
  m4_p$Root_trait <- Root_trait
  Res_em <- rbind(Res_em,m4_em)
  Res_p <- rbind(Res_p,m4_p)
}


#Arrange the pvalues dataframe
Res_p$contrast <- Res_p$contrast %>% gsub(pattern = " ",replacement = "")

#Now the contrast for Mono temperature
c9 <- Res_p$contrast %>% grep(pattern = "Syncom",value = T) %>% unique


####Subset the Mono temperature contrast
Res_p_Mono <- which(Res_p$contrast %in% c9) %>%
  Res_p[.,] %>% droplevels


####Adjust the p-value
Res_p_Mono$Significance <- "NoSignificant"
pval_thres <- 0.1
Res_p_Mono$Significance[which(Res_p_Mono$p.value < pval_thres)] <- "q < 0.1"
Res_p_Mono$Significance <- Res_p_Mono$Significance %>% factor

df_temp_Mono <- Res_p_Mono$contrast %>% gsub(pattern = "-Syncom",replacement = "") %>%
  gsub(pattern = "()", replacement = "") %>%
  as.data.frame

df_temp_Mono

df_temp_Mono$. <- df_temp_Mono$. %>% gsub(pattern = " ",replacement = "")


colnames(df_temp_Mono) <- c("Syncom")

Res_p_Mono <- cbind(Res_p_Mono,df_temp_Mono)

Res_p_Mono$UId <- paste(Res_p_Mono$Genotype,Res_p_Mono$Root_trait,sep = "_")


df_temp_Mono <- data.frame(Syncom = Res_p_Mono$Syncom,
                           p.value =Res_p_Mono$p.value, 
                           estimate = Res_p_Mono$estimate)


df_temp_Mono$Significance <- "NoSignificant"
df_temp_Mono$Significance[which(df_temp_Mono$p.value < pval_thres)] <- "q < 0.1"
df_temp_Mono$Significance <- df_temp_Mono$Significance %>% factor

colnames(df_temp_Mono)='Syncom'

unique(Res_p_Mono$Root_trait)

Res_p_Mono$Syncom = factor(Res_p_Mono$Syncom, c('R289','R245',"NB","P9",'P153','P234',"P129","P118A","P162","R34","R186","P128A.22","P74A","R148A","R345.2","R151","R28","P69","P192","R88","P65"))

ggplot(data = Res_p_Mono, aes(Syncom,Root_trait)) +
  geom_raster(aes(fill = estimate)) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.5,width = 0.85,height = 0.85) + 
  #geom_text(aes(label = Label),color = "black",size = 5)+
  #facet_grid(~Batch,space = "free",scales = "free") +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-3.5,3.5),na.value = "#D9D9D9",oob = squish,name = "Abundance \nz-score") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs. Complete Syncom: ") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.text.x = element_blank()#,
        #legend.text = element_text(size = 5)
  ) +
  scale_y_discrete(expand = c(0,0))
