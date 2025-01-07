##
library(ggplot2)
library(ggpubr)
theme_set(theme_minimal())
library(data.table)
library(ggpmisc)
##
spcs_data<- fread("species.txt")
# summary(lm(CTmax~Tsd, data=spcs_data))
# summary(lm(CTmin~Tsd, data=spcs_data))
cat = as.numeric(unlist(dimnames(table(spcs_data$Tsd))))
data_low<- spcs_data[Tsd==paste(cat[1])]
data_high<- spcs_data[Tsd==paste(cat[11])]
##
a<- ggplot(spcs_data, aes(y = CTmax, x = Tsd))+
  xlab("Long-term variability (°C)")+ ylab("Critical thermal maximum (°C)")+
  coord_cartesian(ylim = c(-20, 40))+
  geom_smooth(method='lm', formula= y~x, color="darkgray", alpha= 0.5)+
  geom_point(stroke=0, alpha= 0.25, aes(color = Tsd))+
  scale_color_gradient(low= "#005CAF", high="#CB1B45")+
  theme(legend.position = "none")

d<- ggplot(spcs_data, aes(y = CTmin, x = Tsd))+
  xlab("Long-term variability (°C)")+ ylab("Critical thermal minimum (°C)")+
  coord_cartesian(ylim = c(-20, 40))+
  geom_smooth(method='lm', formula= y~x, color="darkgray", alpha= 0.5)+
  geom_point(stroke=0, alpha= 0.25, aes(color = Tsd))+
  scale_color_gradient(low= "#005CAF", high="#CB1B45")+
  theme(legend.position = "none")

## Functions for plotting performance curve ####
pcurv<- function(x, x_opt, x_max, sigma, scale){
  if(x< x_opt) temp= (exp(-((x-x_opt)/2/sigma)^2))*scale
  if(x< x_max && x>= x_opt) temp= (1-((x-x_opt)/(x_opt-x_max))^2)*scale
  if(x> x_max) temp= 0
  return(temp)
}

PlotAllLayers<- function(df, col){
  p<- ggplot(data=df, aes(x=df[,1], y=df[,2]))+
    geom_line(alpha= 0.1, color= col)+
    coord_cartesian(xlim = c(-10, 40), ylim = c(0, 3.5))+
    xlab("Temperature (°C)")+ ylab("Performance")
  for(i in 3:length(df[1,])){ 
    p<- p+ geom_line(y=df[,i], alpha= 0.05, color= col)
  }
  return(p)
}
## Making the figure ####
length_iteration<- 500
plot_data<- data.frame()
x_samp<-seq(-10, 40, length.out= 100)
plot_data<- rbind(plot_data,x_samp)
for(i in 1: length_iteration){
  idx<- sample(nrow(data_low),1)
  x_opt= data_low[idx,Topt]
  x_max= data_low[idx,CTmax]
  sigma= data_low[idx,Sigma]
  scale= data_low[idx, Scale]
  temp=rep(NA, length(x_samp))
  for(j in 1: length(x_samp)) temp[j]<- pcurv(x_samp[j],x_opt, x_max, sigma, scale)
  plot_data<- rbind(plot_data, temp)
}

plot_data<- t(plot_data)
plot_data<- as.data.frame(plot_data)
b<- PlotAllLayers(plot_data, col="#005CAF")
#
plot_data<- data.frame()
x_samp<-seq(-10, 40, length.out= 100)
plot_data<- rbind(plot_data,x_samp)
for(i in 1: length_iteration){
  idx<- sample(nrow(data_high),1)
  x_opt= data_high[idx,Topt]
  x_max= data_high[idx,CTmax]
  sigma= data_high[idx,Sigma]
  scale= data_high[idx, Scale]
  temp=rep(NA, length(x_samp))
  for(j in 1: length(x_samp)) temp[j]<- pcurv(x_samp[j],x_opt, x_max, sigma, scale)
  plot_data<- rbind(plot_data, temp)
}

plot_data<- t(plot_data)
plot_data<- as.data.frame(plot_data)
e<- PlotAllLayers(plot_data, col="#CB1B45")
##
xmax<- 50
ymax<- 22.5
#
c<- ggplot(data_low, aes(x=(CTmax-CTmin), y=BRate))+
  coord_cartesian(xlim=c(0, xmax), ylim = c(0, ymax))+
  xlab("Thermal tolerance range (°C)")+ ylab("Species birth rate")+
  geom_point(stroke=0, alpha= 0.5, color= "#005CAF")
#
f<- ggplot(data_high, aes(x=(CTmax-CTmin), y=BRate))+
  coord_cartesian(xlim=c(0, xmax), ylim = c(0, ymax))+
  xlab("Thermal tolerance range (°C)")+ ylab("Species birth rate")+
  geom_point(stroke=0, alpha= 0.5, color= "#CB1B45")
##
fig_comp1<- ggarrange(b,c,nrow=1,labels = c("b","c"))
fig_comp1<- annotate_figure(fig_comp1, top= text_grob("Long-term variability= 0°C (N= 2248)", size=11))
fig_comp2<- ggarrange(e,f,nrow=1,labels = c("e","f"))
fig_comp2<- annotate_figure(fig_comp2, top= text_grob("Long-term variability= 5°C (N= 1366)", size=11))
fig_merge<- ggarrange(a,fig_comp1,d,fig_comp2,nrow=2,ncol=2,widths = c(1,2),labels=c("a","","d",""))

ggsave("LVarExt.pdf",fig_merge, width = 10, height= 6, units="in")
