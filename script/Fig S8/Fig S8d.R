library(ggplot2)

Thicknes = c(68.9, 213, 289.4)
SD = c(3.8, 23.2,43.5)
SR = c('High', 'Medium', 'Low')
Table = cbind(Thicknes, SD)
Table2 = data.frame(cbind(SR, Table))
Table2$Thicknes = as.numeric(Table2$Thicknes)
Table2$SD = as.numeric(Table2$SD)

Table4 = data.frame(cbind(SD, SD))
Table5 = data.frame(cbind(SR, Table4))
paleta_alive <- c("#C200FF",'#FF0000','#8B7500')

Table2$SR = factor(Table2$SR, rev(c('High', 'Medium', 'Low')))

ggplot(data=Table2, aes(x=SR, y=Thicknes, fill=SR)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values = paleta_alive) +
  theme(legend.position="none")+
  theme_few()+
  labs(y='Thicknes (µm)', x='Synthetic Root')+
  geom_pointrange(data = Table5, aes(ymin = Table2$Thicknes+SD,
                  ymax = Table2$Thicknes-SD.1))
