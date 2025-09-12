library(ggplot2)

setwd("G:/My Drive/labs/Nottingham/Duckweed/Exp1/")

Growth2 <- read.table("Data for curves.txt", header = TRUE, row.names = 1)

#Growth curves
LT9243 <- ggplot(Growth2[Growth2$Species=='LT9243',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LT9243',  y='cm2')

LP8539 <- ggplot(Growth2[Growth2$Species=='LP8539',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LP8539',  y='cm2') 

LA7339 <- ggplot(Growth2[Growth2$Species=='LA7339',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LA7339',  y='cm2') 

LM7200 <- ggplot(Growth2[Growth2$Species=='LM7200',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LM7200',  y='cm2') 

LY5280 <- ggplot(Growth2[Growth2$Species=='LY5280',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LY5280',  y='cm2') 

LV7650 <- ggplot(Growth2[Growth2$Species=='LV7650',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LV7650',  y='cm2') 


LM8389 <- ggplot(Growth2[Growth2$Species=='LM8389',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LM8389',  y='cm2') 

LM7123 <- ggplot(Growth2[Growth2$Species=='LM7123',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LM7123',  y='cm2') 

LJ9250 <- ggplot(Growth2[Growth2$Species=='LJ9250',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LJ9250',  y='cm2') 

LT9109 <- ggplot(Growth2[Growth2$Species=='LT9109',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LT9109',  y='cm2') 

LP0049 <- ggplot(Growth2[Growth2$Species=='LP0049',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LP0049',  y='cm2') 

LP7760 <- ggplot(Growth2[Growth2$Species=='LP7760',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='LP7760',  y='cm2') 

SP9192 <- ggplot(Growth2[Growth2$Species=='SP9192',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='SP9192',  y='cm2') 

SP7498 <- ggplot(Growth2[Growth2$Species=='SP7498',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='SP7498',  y='cm2') 

SP5543 <- ggplot(Growth2[Growth2$Species=='SP5543',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='SP5543',  y='cm2') 

SP9509 <- ggplot(Growth2[Growth2$Species=='SP9509',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='SP9509',  y='cm2') 

Pistia <- ggplot(Growth2[Growth2$Species=='PS',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='PS',  y='cm2') 

SI7820 <- ggplot(Growth2[Growth2$Species=='SP7820',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='SI7820',  y='cm2') 

SI9227 <- ggplot(Growth2[Growth2$Species=='SI9227',], aes(x=Days, y=Measure, fill=Treatment)) + 
  geom_point(aes(col=Treatment)) +
  theme_few()+ 
  theme(legend.position="none") +   # draw points
  geom_smooth(method="loess", se=T, level = 0.8) + 
  labs(title='SI9227',  y='cm2')


ggarrange(Pistia,LY5280,SI9227,SP7498,SP9509,
          SI7820,SP5543,SP9192,LM7200,LT9243,
          LP7760,LV7650,LJ9250,LP0049,LM8389,LT9109,
          LA7339,LM7123,LP8539, 
          ncol = 5, nrow=4, title = 'Growth curves',
          font.label = list(size = 6, face = "bold", color ="red"))
