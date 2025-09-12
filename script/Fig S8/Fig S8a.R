###Heatmap
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


setwd("G:/My Drive/labs/Nottingham/Duckweed/thermodynamics/")

Curves <- read.table("Themodynamics2.txt", header = TRUE, row.names = 1)
Curves[2,]

Curves2 = Curves[Curves$ID!=58,]
Curves3 = Curves2[Curves2$ID== 23,]
Curves4 = Curves3[Curves3$ID!= 26,]
Curves5 = Curves4[Curves4$ID!= 27,]

Curves3[Curves3$ID==25,]
Curves3[Curves3$ID==26,]

Curves.mean <- Curves5 %>% group_by(Species, Treatment,Root_complexity,Temperature)%>%
  summarise_all(mean)
Curves.mean2 = data.frame(Curves.mean[Curves.mean$Temperature > 25,])
Curves.mean3 = data.frame(Curves.mean2[Curves.mean2$Temperature < 85,])

LA7339 = ggplot(Curves.mean3[Curves.mean3$Species=='LA7339',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="LA7339") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

LJ9250 = ggplot(Curves.mean3[Curves.mean3$Species=='LJ9250',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="LJ9250") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

LP0049 = ggplot(Curves.mean3[Curves.mean3$Species=='LP0049',], aes(Temperature, Value, fill=ID)) + #, shape=ASV 
  labs(title="LP0049") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')
ylim(-0.1,1)

LP8539 = ggplot(Curves.mean3[Curves.mean3$Species=='LP8539',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="LP8539") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

LT9243 = ggplot(Curves.mean3[Curves.mean3$Species=='LT9243',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="LT9243") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

LY5280 = ggplot(Curves.mean3[Curves.mean3$Species=='LY5280',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="LY5280") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

SI7820 = ggplot(Curves.mean3[Curves.mean3$Species=='SI7820',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="SI7820") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

SI9227 = ggplot(Curves.mean3[Curves.mean3$Species=='SI9227',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="SI9227") +
  theme_bw()+
  theme(legend.position="none") +
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

SP7498 = ggplot(Curves.mean3[Curves.mean3$Species=='SP7498',], aes(Temperature, Value, fill=Treatment)) + #, shape=ASV 
  labs(title="SP7498") +
  theme_bw()+
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)


SP9509 = ggplot(Curves.mean3[Curves.mean3$Species=='SP9509',], aes(Temperature, Value, fill=ID)) + #, shape=ASV 
  labs(title="SP9509") +
  theme_bw()+ 
  theme(legend.position="none")+
  geom_jitter(aes(col=Treatment), size = 1.5)+ 
  labs(x='Temperature',y='mW/g')+
  ylim(-0.1,1)

ggarrange(LA7339,LP0049,LP8539,LJ9250,LT9243,LY5280,SI7820,SI9227,SP7498,SP9509,
          ncol =5, nrow=2)

