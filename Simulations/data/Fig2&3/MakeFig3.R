##
library(ggplot2)
library(ggpubr)
theme_set(theme_minimal())
library(data.table)
library(ggpmisc)
## Functions ####
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
## Median calculation for the entire dataset
# specs<- fread("spcs_summary.txt")
# median(specs$CTmin)
# [1] 3.185495
# median(specs$CTmax)
# [1] 25.89035
## Global options ####
axis_size<- 12
text_size<- 11
## Mean temperature section ####
spcs_data<- fread("mean_svar2.5_lvar0.0/species.txt")
hVol<- fread("mean_svar2.5_lvar0.0/HyperVolume.txt")
##
cat = as.numeric(unlist(dimnames(table(spcs_data$Tmean))))
data_low<- spcs_data[Tmean==paste(cat[1])]
data_high<- spcs_data[Tmean==paste(cat[11])]
##
contours<- make_contours(x=data_low$CTmin, y= data_low$CTmax)
a<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_hline(yintercept = 25.89035, linetype="dashed")+
  geom_vline(xintercept = 3.185495, linetype="dashed")+
  geom_point(data= data_low, aes(y=CTmax, x=CTmin), color="#005CAF", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", linewidth=1)

contours<- make_contours(x=data_high$CTmin, y= data_high$CTmax)
c<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_hline(yintercept = 25.89035, linetype="dashed")+
  geom_vline(xintercept = 3.185495, linetype="dashed")+
  geom_point(data= data_high, aes(y=CTmax, x=CTmin), color="#CB1B45", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", linewidth=1)

## Histogram of TTrange
xmax<- 50
ymax<- 450
binwd<- 3.5
#
data= data_low
dat<- data.frame(x= data$CTmin, y= data$CTmax)
d_out<- dat ######
##
N_all<- nrow(d_out)
N_hot<- nrow(d_out %>% dplyr::filter(x>=3.185495, y>=25.89035))
N_spc<- nrow(d_out %>% dplyr::filter(x>=3.185495, y<25.89035))
N_gen<- nrow(d_out %>% dplyr::filter(x<3.185495, y>=25.89035))
N_cod<- nrow(d_out %>% dplyr::filter(x<3.185495, y<25.89035))
##
plot_data<- data.frame(
  Type=c("Generalist", "Hot specialist","Specialist", "Cold specailist"),
  Proportion=c(N_gen/N_all,N_hot/N_all,N_spc/N_all,N_cod/N_all)
)
##
level_order <- c("Generalist", "Hot specialist","Specialist", "Cold specailist")
b<- ggplot(plot_data, aes(x= Type, y= Proportion))+
  geom_bar(stat="identity", width= 0.2, fill= "#005CAF")+
  coord_cartesian(ylim = c(0, 0.75))+
  scale_x_discrete(limits = level_order)+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(size= axis_size), 
        axis.text.x= element_text(size= text_size))
#
data= data_high
dat<- data.frame(x= data$CTmin, y= data$CTmax)
d_out<- dat ######
##
N_all<- nrow(d_out)
N_hot<- nrow(d_out %>% dplyr::filter(x>=3.185495, y>=25.89035))
N_spc<- nrow(d_out %>% dplyr::filter(x>=3.185495, y<25.89035))
N_gen<- nrow(d_out %>% dplyr::filter(x<3.185495, y>=25.89035))
N_cod<- nrow(d_out %>% dplyr::filter(x<3.185495, y<25.89035))
##
plot_data<- data.frame(
  Type=c("Generalist", "Hot specialist","Specialist", "Cold specailist"),
  Proportion=c(N_gen/N_all,N_hot/ N_all,N_spc/N_all,N_cod/N_all)
)
##
d<- ggplot(plot_data, aes(x= Type, y= Proportion))+
  geom_bar(stat="identity", width= 0.2, fill= "#CB1B45")+
  coord_cartesian(ylim = c(0, 0.75))+
  scale_x_discrete(limits = level_order)+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(size= axis_size), 
        axis.text.x= element_text(size= text_size))

## Short-term variation section ####
spcs_data<- fread("svar_mean18.0_lvar0.0/species.txt")
hVol<- fread("svar_mean18.0_lvar0.0/HyperVolume.txt")
##
cat = as.numeric(unlist(dimnames(table(spcs_data$Tsd_short))))
data_low<- spcs_data[Tsd_short==paste(cat[1])]
data_high<- spcs_data[Tsd_short==paste(cat[11])]
##
contours<- make_contours(x=data_low$CTmin, y= data_low$CTmax)
e<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_hline(yintercept = 25.89035, linetype="dashed")+
  geom_vline(xintercept = 3.185495, linetype="dashed")+
  geom_point(data= data_low, aes(y=CTmax, x=CTmin), color="#005CAF", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", linewidth=1)

contours<- make_contours(x=data_high$CTmin, y= data_high$CTmax)
g<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_hline(yintercept = 25.89035, linetype="dashed")+
  geom_vline(xintercept = 3.185495, linetype="dashed")+
  geom_point(data= data_high, aes(y=CTmax, x=CTmin), color="#CB1B45", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", linewidth=1)
#
data= data_low
dat<- data.frame(x= data$CTmin, y= data$CTmax)
d_out<- dat ######
##
N_all<- nrow(d_out)
N_hot<- nrow(d_out %>% dplyr::filter(x>=3.185495, y>=25.89035))
N_spc<- nrow(d_out %>% dplyr::filter(x>=3.185495, y<25.89035))
N_gen<- nrow(d_out %>% dplyr::filter(x<3.185495, y>=25.89035))
N_cod<- nrow(d_out %>% dplyr::filter(x<3.185495, y<25.89035))
##
plot_data<- data.frame(
  Type=c("Generalist", "Hot specialist","Specialist", "Cold specailist"),
  Proportion=c(N_gen/N_all,N_hot/ N_all,N_spc/N_all,N_cod/N_all)
)
##
f<- ggplot(plot_data, aes(x= Type, y= Proportion))+
  geom_bar(stat="identity", width= 0.2, fill= "#005CAF")+
  coord_cartesian(ylim = c(0, 0.75))+
  scale_x_discrete(limits = level_order)+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(size= axis_size), 
        axis.text.x= element_text(size= text_size))
##
data= data_high
dat<- data.frame(x= data$CTmin, y= data$CTmax)
d_out<- dat ######
##
N_all<- nrow(d_out)
N_hot<- nrow(d_out %>% dplyr::filter(x>=3.185495, y>=25.89035))
N_spc<- nrow(d_out %>% dplyr::filter(x>=3.185495, y<25.89035))
N_gen<- nrow(d_out %>% dplyr::filter(x<3.185495, y>=25.89035))
N_cod<- nrow(d_out %>% dplyr::filter(x<3.185495, y<25.89035))
##
plot_data<- data.frame(
  Type=c("Generalist", "Hot specialist","Specialist", "Cold specailist"),
  Proportion=c(N_gen/N_all,N_hot/ N_all,N_spc/N_all,N_cod/N_all)
)
##
h<- ggplot(plot_data, aes(x= Type, y= Proportion))+
  geom_bar(stat="identity", width= 0.2, fill= "#CB1B45")+
  coord_cartesian(ylim = c(0, 0.75))+
  scale_x_discrete(limits = level_order)+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(size= axis_size), 
        axis.text.x= element_text(size= text_size))

## Long-term variation section ####
spcs_data<- fread("lvar_mean18.0_svar2.5/species.txt")
hVol<- fread("lvar_mean18.0_svar2.5/HyperVolume.txt")
##
cat = as.numeric(unlist(dimnames(table(spcs_data$Tsd))))
data_low<- spcs_data[Tsd==paste(cat[1])]
data_high<- spcs_data[Tsd==paste(cat[11])]
##
contours<- make_contours(x=data_low$CTmin, y= data_low$CTmax)
i<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_hline(yintercept = 25.89035, linetype="dashed")+
  geom_vline(xintercept = 3.185495, linetype="dashed")+
  geom_point(data= data_low, aes(y=CTmax, x=CTmin), color="#005CAF", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", linewidth=1)
##
contours<- make_contours(x=data_high$CTmin, y= data_high$CTmax)
k<- ggplot(NULL)+
  coord_cartesian(xlim=c(-20, 25), ylim = c(15, 45))+
  xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
  theme(axis.title = element_text(size=axis_size))+
  geom_hline(yintercept = 25.89035, linetype="dashed")+
  geom_vline(xintercept = 3.185495, linetype="dashed")+
  geom_point(data= data_high, aes(y=CTmax, x=CTmin), color="#CB1B45", alpha= 0.5, stroke= 0)+
  geom_path(aes(x, y), data=contours$c_50, color="gray25", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_95, color="gray50", linewidth=1) +
  geom_path(aes(x, y), data=contours$c_99, color="gray75", linewidth=1)
#
data= data_low
dat<- data.frame(x= data$CTmin, y= data$CTmax)
d_out<- dat ######
##
N_all<- nrow(d_out)
N_hot<- nrow(d_out %>% dplyr::filter(x>=3.185495, y>=25.89035))
N_spc<- nrow(d_out %>% dplyr::filter(x>=3.185495, y<25.89035))
N_gen<- nrow(d_out %>% dplyr::filter(x<3.185495, y>=25.89035))
N_cod<- nrow(d_out %>% dplyr::filter(x<3.185495, y<25.89035))
##
plot_data<- data.frame(
  Type=c("Generalist", "Hot specialist","Specialist", "Cold specailist"),
  Proportion=c(N_gen/N_all,N_hot/ N_all,N_spc/N_all,N_cod/N_all)
)
##
j<- ggplot(plot_data, aes(x= Type, y= Proportion))+
  geom_bar(stat="identity", width= 0.2, fill= "#005CAF")+
  coord_cartesian(ylim = c(0, 0.75))+
  scale_x_discrete(limits = level_order)+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(size= axis_size), 
        axis.text.x= element_text(size= text_size))
##
data= data_high
dat<- data.frame(x= data$CTmin, y= data$CTmax)
d_out<- dat ######
##
N_all<- nrow(d_out)
N_hot<- nrow(d_out %>% dplyr::filter(x>=3.185495, y>=25.89035))
N_spc<- nrow(d_out %>% dplyr::filter(x>=3.185495, y<25.89035))
N_gen<- nrow(d_out %>% dplyr::filter(x<3.185495, y>=25.89035))
N_cod<- nrow(d_out %>% dplyr::filter(x<3.185495, y<25.89035))
##
plot_data<- data.frame(
  Type=c("Generalist", "Hot specialist","Specialist", "Cold specailist"),
  Proportion=c(N_gen/N_all,N_hot/ N_all,N_spc/N_all,N_cod/N_all)
)
##
l<- ggplot(plot_data, aes(x= Type, y= Proportion))+
  geom_bar(stat="identity", width= 0.2, fill= "#CB1B45")+
  coord_cartesian(ylim = c(0, 0.75))+
  scale_x_discrete(limits = level_order)+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(size= axis_size), 
        axis.text.x= element_text(size= text_size))

#
# Merge figures ####
title_size<- 16
fig_comp1<- ggarrange(a,b,nrow=1,labels = c("a","b"))
fig_comp1<- annotate_figure(fig_comp1, top= text_grob("Low ambient temperature (N= 1776)", size=title_size))
fig_comp2<- ggarrange(c,d,nrow=1,labels = c("c","d"))
fig_comp2<- annotate_figure(fig_comp2, top= text_grob("High ambient temperature (N= 2198)", size=title_size))
fig_row1<- ggarrange(fig_comp1,fig_comp2,nrow=1)
#
fig_comp1<- ggarrange(e,f,nrow=1,labels = c("e","f"))
fig_comp1<- annotate_figure(fig_comp1, top= text_grob("Low short-term variability (N= 2233)", size=title_size))
fig_comp2<- ggarrange(g,h,nrow=1,labels = c("g","h"))
fig_comp2<- annotate_figure(fig_comp2, top= text_grob("High short-term variability (N= 2368)", size=title_size))
fig_row2<- ggarrange(fig_comp1,fig_comp2,nrow=1)
#
fig_comp1<- ggarrange(i,j,nrow=1,labels = c("i","j"))
fig_comp1<- annotate_figure(fig_comp1, top= text_grob("Low long-term variability (N= 2248)", size=title_size))
fig_comp2<- ggarrange(k,l,nrow=1,labels = c("k","l"))
fig_comp2<- annotate_figure(fig_comp2, top= text_grob("High long-term variability (N= 1366)", size=title_size))
fig_row3<- ggarrange(fig_comp1,fig_comp2,nrow=1)
#
fig_merge<- ggarrange(fig_row1,fig_row2,fig_row3, ncol=1)
##
ggsave("Fig3.pdf",fig_merge, width = 16.5, height= 12.5, units="in")
