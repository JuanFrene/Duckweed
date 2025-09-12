library(ggplot2)
library(scales)
library(emmeans)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Root traits/Exp1/")

TablaCompleta = read.table("Root traits.txt", header = TRUE, row.names = 1)

TablaCompleta$Microbiome = factor(TablaCompleta$Microbiome, c('NB', 'B+'))

TablaCompleta$Species = factor(TablaCompleta$Species, c('PS','LY5280','SI9227','SP7498','SP9509',
                                                        'SI7820','SP5543','SP9192','LM7200','LT9243',
                                                        'LP7760','LV7650','LJ9250','LP0049','LM8389','LT9109',
                                                        'LA7339','LM7123','LP8539'))

#Do a model per Root complexity
Res_em <- NULL
Res_p <- NULL
for(specie in TablaCompleta$Species %>% unique){
  melted_sub <- TablaCompleta %>% subset(Species  == specie) %>% droplevels
  m1 <- lm(data = melted_sub,
           formula = Stile_diameter ~ Microbiome)
  m1_res <- emmeans(m1,pairwise~Microbiome,ref =1, adjust = "none")
  m1_em <- m1_res$emmeans %>% as.data.frame
  m1_p <- m1_res$contrasts %>% as.data.frame
  m1_em$Root_Trait <- specie
  m1_p$Root_Trait <- specie
  Res_em <- rbind(Res_em,m1_em)
  Res_p <- rbind(Res_p,m1_p)
}
Res_p

S_C_NB_Mean = aggregate(CCL2__n_cell  ~Species,TablaCompleta[TablaCompleta$Microbiome=='NB',],mean)
S_C_B_Mean = aggregate(CCL2__n_cell  ~Species,TablaCompleta[TablaCompleta$Microbiome=='B+',],mean)
Tabla=cbind(S_C_NB_Mean,S_C_B_Mean[,2])
colnames(Tabla)=c('Species','NB','B')
Tabla$NB = as.numeric(Tabla$NB)
Tabla$B = as.numeric(Tabla$B)

# prep data
left_label <- paste(Tabla$Species, round(Tabla$NB),sep=", ")
right_label <- paste(Tabla$Species, round(Tabla$B),sep=", ")
Tabla$class <- ifelse((Tabla$NB - Tabla$B) <= 0,  "green","red")

# Plot
# Plot
p <- ggplot(Tabla) + geom_segment(aes(x=1, xend=2, y=NB, yend=B), size=.75, show.legend=F) + 
  geom_vline(xintercept=1, linetype="dashed", linewidth=.05) + 
  geom_vline(xintercept=2, linetype="dashed", linewidth=.05) +
  scale_color_manual(labels = c("Up", "Down"), 
                     values = c("green"="#00ba38", "red"="#f8766d")) +  # color of lines
  labs(x="", y="CCL2 N cell ")  # Axis labels
  xlim(.5, 2.5) + ylim(5,(1.1*(max(Tabla$NB, Tabla$B))))  # X and Y axis limits

# Add texts
p <- p + geom_text(label=left_label, y=Tabla$NB, x=rep(1, NROW(Tabla)), hjust=1.1, size=3.5)
p <- p + geom_text(label=right_label, y=Tabla$B, x=rep(2, NROW(Tabla)), hjust=-0.1, size=3.5)
p <- p + geom_text(label="Sterile", x=1, y=1.1*(max(Tabla$NB, Tabla$B)), hjust=1.2, size=5)  # title
p <- p + geom_text(label="Natural", x=2, y=1.1*(max(Tabla$NB, Tabla$B)), hjust=-0.1, size=5)  # title

# Minify theme
p + theme(panel.background = element_blank(), 
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"))

