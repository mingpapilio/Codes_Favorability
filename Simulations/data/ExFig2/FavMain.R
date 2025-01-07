##
library(ggplot2)
library(ggpubr)
theme_set(theme_minimal())
library(data.table)
library(ggpmisc)

##
spcs_data<- fread("species.txt")
hVol<- fread("HyperVolume.txt")
##
a<- ggplot(hVol, aes(y=CTmax_CTmin_hv, x= Tmean))+
  xlab("Mean ambient temperature (°C)")+ ylab("Hypervolume of functional traits")+
  coord_cartesian(ylim = c(9, 17))+
  geom_smooth(method='lm', formula= y~x, color="darkgray", alpha= 0.5)+
  geom_point(stroke=0, size= 3, aes(color = Tmean))+
  scale_color_gradient(low="#005CAF", high="#CB1B45")+
  theme(legend.position="none")+
  stat_poly_eq(aes(label= paste(
    after_stat(eq.label),
    after_stat(p.value.label), 
    after_stat(rr.label), sep = "*\", \"*")))

##
cat = as.numeric(unlist(dimnames(table(spcs_data$Tmean))))
data_low<- spcs_data[Tmean==paste(cat[1])]
data_high<- spcs_data[Tmean==paste(cat[11])]

make_contours<- function(x, y){
  d<- data.frame(x= x, y= y)
  kd<- ks::kde(d, compute.cont=TRUE)
  contour_99<- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["1%"])[[1]])
  contour_99<- data.frame(contour_99)
  contour_95<- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["5%"])[[1]])
  contour_95<- data.frame(contour_95)
  contour_50<- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["50%"])[[1]])
  contour_50<- data.frame(contour_50)
  out<- list(c_50=contour_50,c_95=contour_95,c_99=contour_99)
  return(out)
}
##
contours<- make_contours(x=data_low$Tmin, y= data_low$Tmax)
b<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=8))+
  geom_point(data= data_low, aes(y=Tmax, x=Tmin), color="#005CAF", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", size=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", size=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", size=1)

contours<- make_contours(x=data_high$Tmin, y= data_high$Tmax)
d<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=8))+
  geom_point(data= data_high, aes(y=Tmax, x=Tmin), color="#CB1B45", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", size=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", size=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", size=1)

## Histogram of TTrange
xmax<- 50
ymax<- 450
binwd<- 3.5
#
c<- ggplot(data_low,aes(x=(Tmax-Tmin)))+
  coord_cartesian(xlim=c(0, xmax), ylim = c(0, ymax))+
  xlab("Thermal tolerance range (°C)")+ ylab("Number of species")+
  theme(axis.title = element_text(size=8))+
  geom_histogram(binwidth= binwd, na.rm = T, fill="#005CAF", color= "white")
#
e<- ggplot(data_high,aes(x=(Tmax-Tmin)))+
  coord_cartesian(xlim=c(0, xmax), ylim = c(0, ymax))+
  xlab("Thermal tolerance range (°C)")+ ylab("Number of species")+
  theme(axis.title = element_text(size=8))+
  geom_histogram(binwidth= binwd, na.rm = T, fill= "#CB1B45", color= "white")
#
fig_comp1<- ggarrange(b,c,nrow=1,labels = c("b","c"))
fig_comp1<- annotate_figure(fig_comp1, top= text_grob("Mean ambient temperature= 14°C (N= 2253)", size=11))
fig_comp2<- ggarrange(d,e,nrow=1,labels = c("d","e"))
fig_comp2<- annotate_figure(fig_comp2, top= text_grob("Mean ambient temperature= 22°C (N= 2232)", size=11))
fig_comp3<- ggarrange(fig_comp1, fig_comp2, nrow=2)
fig_merge<- ggarrange(a,fig_comp3, nrow=1,labels=c("a",""))
##
ggsave("FavMain.pdf",fig_merge, width = 10, height= 5, units="in")
 