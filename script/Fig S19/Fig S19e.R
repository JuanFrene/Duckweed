#### DAPI cuantification

DAPI <- read.table("DAPI cuantification.txt", header = TRUE, row.names = 1)

DAPI$Root = factor(DAPI$Root, c( 'High', 'Medium','Low'))
DAPI$Species = factor(DAPI$Species, c('LY9205', 'SI9227', 'SP7498', 'SP9509', 'SI7820', 'LT9243', 'LJ9250', 'LP0049', 'LA7339', 'LP8539'))

paleta_alive <- c("#C200FF",'#FF0000','#8B7500','#00008B',"#FFB919",'#FF7F50',"#00CC7A")

DAPI$Complexity = as.numeric(DAPI$Complexity)
theme_set(theme_bw())  # pre-set the bw theme.
ggplot(DAPI, aes(Complexity, Abundance)) +  
  geom_jitter(aes(col=Species, shape=Root), size = 2.5) + 
  geom_smooth(method="lm", se=T, col = 'black')+
  facet_grid(Treatment~.,scales = "fixed",space = "fixed")+
  theme_few()+
  labs(y='Abundance', x= 'Root complexity')
geom_text(aes(label = "p = 2.44e-05",
              x = 7.35, y = 105, size = 4, fill = "black"))+#R2=0.065
  geom_text(aes(label = "R  = 0.1841",
                x = 7.3, y = 110, size = 4, fill = "black"))


summary(lm(Complexity~Abundance, data = DAPI))

