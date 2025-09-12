#####Load packages
library(Rmisc)
library(ggrepel)
library(scales)
library(ggtree)
library(harrietr)
library(emmeans)

######Set work directory
setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp Metabolites addition/")


######Open data file
Data <- read.table("Root trait LC50.txt", header = TRUE, row.names = 1)
head(Data)
ID=rownames(Data)
Data= cbind(ID,Data)
###Transform in z-score matrix 
#A z-score is a measure of how many standard deviations 
#(how far my observation is from the mean) below or above
#the populatio mean a raw score is.
Tab_f <- Data[,-(1:4)] %>% scale
nrow(Tab_f)

### Now we need to summarize 
melted <- Tab_f %>% melt
colnames(melted) <- c("ID","parameter","value")

###melted data and metadata
melted <- merge(Data[,1:4],melted,  by='ID')
melted $ Treat <- paste (melted$Treatment, melted$Dosis, sep = "_")


#Do a model per ion
Res_em <- NULL
Res_p <- NULL

for(ion in melted$parameter %>% unique){
  melted_sub <- melted %>% subset(parameter  == ion) %>% droplevels
  m4 <- lm(data = melted_sub,
           formula = value ~ Treat)
  m4_res <- emmeans(m4,pairwise~Treat,ref =4, adjust = "none")
  m4_em <- m4_res$emmeans %>% as.data.frame
  m4_p <- m4_res$contrasts %>% as.data.frame
  m4_em$Ion <- ion
  m4_p$Ion <- ion
  Res_em <- rbind(Res_em,m4_em)
  Res_p <- rbind(Res_p,m4_p)
}



#Arrange the pvalues dataframe
Res_p$contrast <- Res_p$contrast %>% gsub(pattern = " ",replacement = "")

#Now the contrast for Room temperature
c9 <- Res_p$contrast %>% grep(pattern = "Syncom_Control",value = T) %>% unique


####Subset the Room temperature contrast
Res_p_com <- which(Res_p$contrast %in% c9) %>%
  Res_p[.,] %>% droplevels

nrow(Res_p_com)

Res_p_com$Significance <- "NoSignificant"
pval_thres <- 0.1
Res_p_com$Significance[which(Res_p_com$p.value < pval_thres)] <- "q < 0.1"
Res_p_com$Significance <-Res_p_com$Significance %>% factor

Res_p_com$contrast <- Res_p_com$contrast %>% gsub(pattern = "-Syncom_Control",replacement = "")


Compound=c('C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5'
           ,'C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4'
           ,'C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4',
           'C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB','C4','C4','C4','C5','C5','C5','NB')
Res_p_com2=cbind(Res_p_com,Compound)

Res_p_com2$contrast = factor(Res_p_com2$contrast, c("NB_NB","C4_25nM","C4_50mM","C4_75mM","C5_25nM","C5_50mM","C5_75mM"))
Res_p_com2$Ion= factor(Res_p_com2$Ion, c('CCL','complexity','Steele_cell_n',	'Endodermis_cell_n',	'CCL1_cell_n',	'CCL2_cell_n',	'CCL3_n_cell',	'CCL4_n_cell',	'Exodermis_n_cell','Epidermis_n_cell','Aerenchyma_mean_area','Aerenchyma_n','Aerenchyma_total_area','Steele','Endodermis','CCL1','CCL2','CCL3','CCL4','Exodermis','Epidermis','Max_diameter'))
Res_p_com2$Compound = factor(Res_p_com2$Compound, c("NB","C4","C5"))

ggplot(data = Res_p_com2[Res_p_com2$Ion!='complexity',], aes(contrast,Ion)) +
  geom_raster(aes(fill = estimate)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.5,width = 0.95,height = 0.95) + 
  #geom_text(aes(label = Label),color = "black",size = 5)+
  facet_grid(~Compound,space = "free",scales = "free") +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-3.5,3.5),na.value = "#D9D9D9",oob = squish,name = "Abundance \nz-score") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs. Control ") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.position = "right",
        strip.text.y = element_blank(),
        strip.text.x = element_blank()#,
  ) +
  scale_y_discrete(expand = c(0,0))


