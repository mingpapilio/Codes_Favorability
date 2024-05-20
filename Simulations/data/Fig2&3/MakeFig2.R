##
library(ggplot2)
library(ggpubr)
theme_set(theme_minimal())
library(data.table)
library(ggpmisc)
## Functions
make_contour95<- function(x, y){
  dat<- data.frame(x= x, y= y)
  kd<- ks::kde(dat, compute.cont=TRUE)
  contour_95<- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["5%"])[[1]])
  contour_95<- data.frame(contour_95)
  return(contour_95)
}
## Global options ####
ymax<- 21
ymin<- 11
axis_size<- 12
## Mean temperature section ####
spcs_data<- fread("mean_svar2.5_lvar0.0/species.txt")
hVol<- fread("mean_svar2.5_lvar0.0/HyperVolume.txt")
##
a<- ggplot(hVol, aes(y=CTmax_CTmin_hv, x= Tmean))+
  xlab("Mean ambient temperature (°C)")+ ylab("Hypervolume of functional traits")+
  coord_cartesian(ylim = c(ymin, ymax))+
  geom_smooth(method='lm', formula= y~x, color="darkgray", alpha= 0.5)+
  geom_point(stroke=0, size= 3, aes(color = Tmean))+
  scale_color_gradient(low="#005CAF", high="#CB1B45")+
  theme(legend.position="none", axis.title = element_text(size=axis_size))+
  stat_poly_eq(aes(label= paste(
    after_stat(eq.label),
    after_stat(p.value.label), 
    after_stat(rr.label), sep = "*\", \"*")))
##
cat = as.numeric(unlist(dimnames(table(spcs_data$Tmean))))
#
contours<- list()
for (i in c(1,6,11)){
  focal<- spcs_data[Tmean==paste(cat[i])]
  contr<- make_contour95(x= focal$CTmin, y=focal$CTmax)
  contours<- append(contours, list(contr))
}
#
b<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_path(aes(x, y), data=contours[[1]], color="#005CAF", linewidth= 1, alpha= 0.5)+
  geom_path(aes(x, y), data=contours[[2]], color="#653B79", linewidth= 1, alpha= 0.5)+
  geom_path(aes(x, y), data=contours[[3]], color="#CB1B45", linewidth= 1, alpha= 0.5)
## Short-term variation section ####
spcs_data<- fread("svar_mean18.0_lvar0.0/species.txt")
hVol<- fread("svar_mean18.0_lvar0.0/HyperVolume.txt")
##
c<- ggplot(hVol, aes(y=CTmax_CTmin_hv, x= Tsd_short))+
  xlab("Short-term variability (°C)")+ ylab("Hypervolume of functional traits")+
  coord_cartesian(ylim = c(ymin, ymax))+
  geom_smooth(method='lm', formula= y~x, color="darkgray", alpha= 0.5)+
  geom_point(stroke=0, size= 3, aes(color = Tsd_short))+
  scale_color_gradient(low="#005CAF", high="#CB1B45")+
  theme(legend.position="none", axis.title = element_text(size=axis_size))+
  stat_poly_eq(aes(label= paste(
    after_stat(eq.label),
    after_stat(p.value.label), 
    after_stat(rr.label), sep = "*\", \"*")))

##
cat = as.numeric(unlist(dimnames(table(spcs_data$Tsd_short))))
contours<- list()
for (i in c(1,6,11)){
  focal<- spcs_data[Tsd_short==paste(cat[i])]
  contr<- make_contour95(x= focal$CTmin, y=focal$CTmax)
  contours<- append(contours, list(contr))
}
#
d<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_path(aes(x, y), data=contours[[1]], color="#005CAF", linewidth= 1, alpha= 0.5)+
  geom_path(aes(x, y), data=contours[[2]], color="#653B79", linewidth= 1, alpha= 0.5)+
  geom_path(aes(x, y), data=contours[[3]], color="#CB1B45", linewidth= 1, alpha= 0.5)
## Long-term variation section ####
spcs_data<- fread("lvar_mean18.0_svar2.5/species.txt")
hVol<- fread("lvar_mean18.0_svar2.5/HyperVolume.txt")
##
e<- ggplot(hVol, aes(y=CTmax_CTmin_hv, x= Tsd))+
  xlab("Long-term variability (°C)")+ ylab("Hypervolume of functional traits")+
  coord_cartesian(ylim = c(ymin, ymax))+
  geom_smooth(method='lm', formula= y~x, color="darkgray", alpha= 0.5)+
  geom_point(stroke=0, size= 3, aes(color = Tsd))+
  scale_color_gradient(low="#005CAF", high="#CB1B45")+
  theme(legend.position="none", axis.title = element_text(size=axis_size))+
  stat_poly_eq(aes(label= paste(
    after_stat(eq.label),
    after_stat(p.value.label), 
    after_stat(rr.label), sep = "*\", \"*")))

##
cat = as.numeric(unlist(dimnames(table(spcs_data$Tsd))))
contours<- list()
for (i in c(1,6,11)){
  focal<- spcs_data[Tsd==paste(cat[i])]
  contr<- make_contour95(x= focal$CTmin, y=focal$CTmax)
  contours<- append(contours, list(contr))
}
#
f<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_path(aes(x, y), data=contours[[1]], color="#005CAF", linewidth= 1, alpha= 0.5)+
  geom_path(aes(x, y), data=contours[[2]], color="#653B79", linewidth= 1, alpha= 0.5)+
  geom_path(aes(x, y), data=contours[[3]], color="#CB1B45", linewidth= 1, alpha= 0.5)
## Merge figure ####
fig_comp1<- ggarrange(a,b,nrow=1,labels = c("a","b"))
fig_comp1<- annotate_figure(fig_comp1, top= text_grob("Mean ambient temperature", size=16))
fig_comp2<- ggarrange(c,d,nrow=1,labels = c("c","d"))
fig_comp2<- annotate_figure(fig_comp2, top= text_grob("Short-term variability", size=16))
fig_comp3<- ggarrange(e,f,nrow=1,labels = c("e","f"))
fig_comp3<- annotate_figure(fig_comp3, top= text_grob("Long-term variability", size=16))
fig_merge<- ggarrange(fig_comp1,fig_comp2,fig_comp3, ncol=1)
##
ggsave("Fig2.pdf",fig_merge, width = 8.5, height= 12.5, units="in")
