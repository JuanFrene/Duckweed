library(ggplot2)
library(cowplot)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Bacrteria Inhibition exp/")

# 1. Metadata
Data <- read.table("Bacteria_inhibition.txt", header = TRUE, row.names = 1)
Slope <- read.table("Bacteria_inhibition_slope.txt", header = TRUE, row.names = 1)

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

Data1 = Data[Data$Experiment == '1',]
Data2 = Data[Data$Experiment != '1',]

log2foldchange2 <- c()
for(specie in Data2$Microbe  %>% unique){
  melted_sub <- Data2 %>% subset(Microbe ==  specie) %>% droplevels
  for(time in melted_sub[melted_sub$Time!='T0',]$Time  %>% unique){
    melted_sub2 <- melted_sub[melted_sub$Time!='T0',] %>% subset(Time ==  time) %>% droplevels #'ASV522'
    for(drug in melted_sub2[melted_sub2$Compound2  =='Compound',]$Compound  %>% unique){
      melted_sub3 <- melted_sub2[melted_sub2$Compound2  =='Compound',] %>% subset(Compound  ==  drug) %>% droplevels
      
      Control = data.frame(t(melted_sub2[melted_sub2$Compound  =='Control',][,10]))
      Compound = data.frame(t(melted_sub3[melted_sub3$Compound2  =='Compound',][,10]))
      
      Control_mean = apply(Control, 1, mean) 
      Compound_mean = apply(Compound, 1, mean) 
      
      C_C_mean <- log2(Compound_mean+1) - log2(Control_mean+1) 
      
      C_C_mean_statistic = t.test(Control,Compound)
      pvalueB_C_C = C_C_mean_statistic$p.value
      C_C_mean_conf = t(data.frame(C_C_mean_statistic$conf.int))
      colnames(C_C_mean_conf) = c('inf','sup')
      
      result = cbind(C_C_mean,pvalueB_C_C)
      specie2 = specie
      row.names(result)= time
      result2 = cbind(specie2,result)
      result3 = cbind(drug,result2)
      result4 = cbind(time,result3)
      result5 = cbind(result4,C_C_mean_conf)
      log2foldchange2 <- rbind(log2foldchange2,result5)
    }}}

log2foldchange_a =data.frame(rbind(log2foldchange1,log2foldchange2))

log2foldchange_a$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange_a$Significance[which(log2foldchange_a$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange_a$Significance <- log2foldchange_a$Significance %>% factor
log2foldchange_a$C_C_mean = as.numeric(log2foldchange_a$C_C_mean)

log2foldchange_a$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange_a$Significance[which(log2foldchange_a$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange_a$Significance <- log2foldchange_a$Significance %>% factor
log2foldchange_a$C_C_mean = as.numeric(log2foldchange_a$C_C_mean)


head(log2foldchange2)

# Plot
log2foldchange = log2foldchange_a[log2foldchange_a$specie2!='Control_Neg',]

#C1
new_df <- log2foldchange[log2foldchange$drug=='C1',] %>% 
  # desc orders from largest to smallest
  arrange(desc(C_C_mean))

new_df2 = data.frame(new_df[new_df$time=='T6',])
order = 1:191
new_df3 = cbind(new_df2, order)
new_df3$specie2 = factor(new_df3$specie2, rev(new_df3$specie2))

C1 = ggplot(new_df3, aes(x=specie2, y=C_C_mean)) + 
  geom_hline(yintercept = 0,size = 0.5,color = "black")+
  geom_point(stat='identity', aes(col=Significance), size=2)  +
  scale_color_manual(values = c('black',"red"))+
  labs(title="Methy-Citrate", y='log2FC', x='Syncom members') + 
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(size = 3),
        legend.position = 'none')

#C2
new_df <- log2foldchange[log2foldchange$drug=='C2',] %>% 
  # desc orders from largest to smallest
  arrange(desc(C_C_mean))

new_df2 = data.frame(new_df[new_df$time=='T6',])
order = 1:191
new_df3 = cbind(new_df2, order)
new_df3$specie2 = factor(new_df3$specie2, rev(new_df3$specie2))

C2 = ggplot(new_df3, aes(x=specie2, y=C_C_mean)) + 
  geom_hline(yintercept = 0,size = 0.5,color = "black")+
  geom_point(stat='identity', aes(col=Significance), size=2)  +
  scale_color_manual(values = c('black',"red"))+
  labs(title="Citruline", y='log2FC', x='Syncom members') + 
  #ylim(-2.5, 2.5) +
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(size = 3),
        legend.position = 'none')

#C3
new_df <- log2foldchange[log2foldchange$drug=='C3',] %>% 
  # desc orders from largest to smallest
  arrange(desc(C_C_mean))

new_df2 = data.frame(new_df[new_df$time=='T6',])
order = 1:191
new_df3 = cbind(new_df2, order)
new_df3$specie2 = factor(new_df3$specie2, rev(new_df3$specie2))

C3 = ggplot(new_df3, aes(x=specie2, y=C_C_mean)) + 
  geom_hline(yintercept = 0,size = 0.5,color = "black")+
  geom_point(stat='identity', aes(col=Significance), size=2)  +
  scale_color_manual(values = c('black',"red"))+
  labs(title="Adrenaline", y='log2FC', x='Syncom members') + 
  #ylim(-2.5, 2.5) +
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(size = 3),
        legend.position = 'none')


#C4
new_df <- log2foldchange[log2foldchange$drug=='C3',] %>% 
  # desc orders from largest to smallest
  arrange(desc(C_C_mean))

new_df2 = data.frame(new_df[new_df$time=='T6',])
order = 1:191
new_df3 = cbind(new_df2, order)
new_df3$specie2 = factor(new_df3$specie2, rev(new_df3$specie2))

C4 = ggplot(new_df3, aes(x=specie2, y=C_C_mean)) + 
  geom_hline(yintercept = 0,size = 0.5,color = "black")+
  geom_point(stat='identity', aes(col=Significance), size=2)  +
  scale_color_manual(values = c('black',"red"))+
  labs(title="3-Methyl-Lysine", y='log2FC', x='Syncom members') + 
  #ylim(-2.5, 2.5) +
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(size = 3),
        legend.position = 'none')

#C5
new_df <- log2foldchange[log2foldchange$drug=='C5',] %>% 
  # desc orders from largest to smallest
  arrange(desc(C_C_mean))

new_df2 = data.frame(new_df[new_df$time=='T6',])
order = 1:191
new_df3 = cbind(new_df2, order)
new_df3$specie2 = factor(new_df3$specie2, rev(new_df3$specie2))

C5 = ggplot(new_df3, aes(x=specie2, y=C_C_mean)) + 
  geom_hline(yintercept = 0,size = 0.5,color = "black")+
  geom_point(stat='identity', aes(col=Significance), size=2)  +
  scale_color_manual(values = c('black',"red"))+
  labs(title="Glutaric acid", y='log2FC', x='Syncom members') + 
  #ylim(-2.5, 2.5) +
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(size = 3),
        legend.position = 'none')

#C6
new_df <- log2foldchange[log2foldchange$drug=='C6',] %>% 
  # desc orders from largest to smallest
  arrange(desc(C_C_mean))

new_df2 = data.frame(new_df[new_df$time=='T6',])
order = 1:191
new_df3 = cbind(new_df2, order)
new_df3$specie2 = factor(new_df3$specie2, rev(new_df3$specie2))

C6 = ggplot(new_df3, aes(x=specie2, y=C_C_mean)) +
  geom_hline(yintercept = 0,size = 0.5,color = "black")+
  geom_point(stat='identity', aes(col=Significance), size=2)  +
  scale_color_manual(values = c('black',"red"))+
  labs(title="Pyrone", y='log2FC', x='Syncom members') + 
  #ylim(-2.5, 2.5) +
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(size = 3),
        legend.position = 'none')


plot_grid(C1,C2,C3,C4,C5,C6, align = 'h', rel_widths = c(1,1,1,1,1,1))
