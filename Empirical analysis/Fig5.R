library(dplyr)
library(tidyr)
library(vegan)

## load data
all_data = df.all = read.csv('./FAV_hv_data.csv',header = TRUE)

## grouping species into high, mid and low STmean, DTR and STR ####
all_data = all_data %>% mutate(Tmean_g = ifelse(Group == 'M1'|Group == 'M2'|Group == 'T1','High', 
                                                ifelse(Group == 'T2'|Group == 'C1'|Group == 'T3', 'Mid', 'Low')))
all_data = all_data %>% mutate(DTR_g = ifelse(Group == 'C1'|Group == 'C3'|Group == 'T1', 'High',
                                              ifelse(Group == 'C2'|Group == 'C4'|Group == 'T2', 'Mid','Low')))
all_data = all_data %>% mutate(STR_g = ifelse(Group == 'C4'|Group == 'C3'|Group == 'C2', 'High',
                                              ifelse(Group == 'C1'|Group == 'T2'|Group == 'T3', 'Mid','Low')))
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

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

## Median calculation for the entire dataset
ctmin_mean = median(all_data$CTmin)
ctmax_mean = median(all_data$CTmax)
## Global options ####
axis_size<- 12
text_size<- 13
#---- Mean temperature section ----
plot_list = list()
i = 1
for(env_type in c('Tmean', 'DTR', 'STR')){
  for(group in c('Low', 'Mid', 'High')){
    if(group == 'High') col = '#CB1B45'
    else if(group == 'Mid') col = '#937DC2'
    else if(group == 'Low') col = '#005CAF'
    if(env_type == 'Tmean') data_low = all_data %>% filter(Tmean_g == group) %>% select(CTmin, CTmax)
    else if(env_type == 'DTR') data_low = all_data %>% filter(DTR_g == group) %>% select(CTmin, CTmax)
    else if(env_type == 'STR') data_low = all_data %>% filter(STR_g == group) %>% select(CTmin, CTmax)
    
    contours<- make_contours(x=data_low$CTmin, y= data_low$CTmax)
    a<- ggplot(NULL)+
      coord_cartesian(xlim=c(-10,20), ylim = c(25,50))+
      xlab("Critical thermal minimum (°C)")+ ylab("Critical thermal maximum (°C)")+
      theme(axis.title = element_text(size=axis_size),
            axis.text = element_text(size=text_size))+
      geom_hline(yintercept = ctmax_mean, linetype="dashed")+
      geom_vline(xintercept = ctmin_mean, linetype="dashed")+
      geom_point(data= data_low, aes(y=CTmax, x=CTmin), color=col, alpha= 0.5, stroke= 0)+
      geom_path(aes(x, y), data=contours$c_50, color="gray25", linewidth=1) +
      geom_path(aes(x, y), data=contours$c_95, color="gray50", linewidth=1) +
      geom_path(aes(x, y), data=contours$c_99, color="gray75", linewidth=1)
    
    ## Histogram of proportion
    #
    data= data_low
    dat<- data.frame(x= data$CTmin, y= data$CTmax)
    d_out=dat
    ##
    N_all<- nrow(d_out)
    N_hot<- nrow(d_out %>% dplyr::filter(x>=ctmin_mean, y>=ctmax_mean))
    N_spc<- nrow(d_out %>% dplyr::filter(x>=ctmin_mean, y<ctmax_mean))
    N_gen<- nrow(d_out %>% dplyr::filter(x<ctmin_mean, y>=ctmax_mean))
    N_cod<- nrow(d_out %>% dplyr::filter(x<ctmin_mean, y<ctmax_mean))
    ##
    plot_data<- data.frame(
      Type=c("Generalist", "Hot-adapted species","Specialist", "Cold-adapted species"),
      Proportion=c(N_gen/N_all,N_hot/N_all,N_spc/N_all,N_cod/N_all)
    )
    ##
    level_order <- c("Generalist", "Hot-adapted species","Specialist", "Cold-adapted species")
    b<- ggplot(plot_data, aes(x= Type, y= Proportion))+
      geom_bar(stat="identity", width= 0.2, fill= col)+
      coord_cartesian(ylim = c(0, 0.75))+
      scale_x_discrete(limits = level_order, 
                       labels=addline_format(c("Generalist", "Hot-adapted species",
                                               "Specialist", "Cold-adapted species")))+
      theme(axis.title.x= element_blank(), 
            axis.title.y= element_text(size= axis_size), 
            axis.text= element_text(size= text_size))
    
      plot_list[[i]] = a
      plot_list[[i+1]] = b
      i = i+2
  }
}

## Merge figures ##
title_size<- 16
main_text_1 = c('Low','Medium','High')
main_text_2 = c('annual mean temperature','diurnal temperature range','seasonal temperature range')

fig_comp = list()
fig_row = list()
i = j = 1
for(m in seq(1,3)){
  for(k in seq(1,3)){
    fig_comp1<- ggarrange(plot_list[[i]],plot_list[[i+1]],nrow=1,labels = c(letters[i],letters[i+1]))
    fig_comp[[j]]<- annotate_figure(fig_comp1, top= text_grob(paste(main_text_1[k], main_text_2[m]), size=title_size))
    i = i+2
    j = j+1
  }
  fig_row[[m]]<- ggarrange(fig_comp[[m*3-2]],fig_comp[[m*3-1]],fig_comp[[m*3]],nrow=1)
}

#
fig_merge<- ggarrange(fig_row[[1]],fig_row[[2]],fig_row[[3]], ncol=1)
##
ggsave("./Fig5_kde.pdf",fig_merge, width = 25, height= 12.5, units="in")

#---- Histogram of TTrange ----
xmax<- 50 # Maximal x value in the plot
ymax<- 100 # Maximal y value in the plot
binwd<- 2 # Width of each bin
axis_size<- 12
text_size<- 13

plot_list2 = list()
i = 1

for(env_type in c('Tmean', 'DTR', 'STR')){
  for(group in c('Low', 'Mid', 'High')){
    if(group == 'High') col = '#CB1B45'
    else if(group == 'Mid') col = '#937DC2'
    else if(group == 'Low') col = '#005CAF'
    if(env_type == 'Tmean') data_low = all_data %>% filter(Tmean_g == group) %>% select(TTrange)
    else if(env_type == 'DTR') data_low = all_data %>% filter(DTR_g == group) %>% select(TTrange)
    else if(env_type == 'STR') data_low = all_data %>% filter(STR_g == group) %>% select(TTrange)
        
    h<- ggplot(data_low, aes(x= TTrange))+
      coord_cartesian(xlim=c(20, 50), ylim = c(0, ymax))+
      ggtitle(paste(group,env_type))+
      xlab("TTrange (°C)")+ ylab("")+
      geom_histogram(binwidth= binwd, na.rm = T, fill=col, color='grey80')+
      theme(axis.title.y= element_text(size= axis_size), 
            axis.text= element_text(size= text_size))
    h
    
    plot_list2[[i]] = h
    i=i+1
  }}

##
library(gridExtra)
hist_merge = grid.arrange(grobs = plot_list2, ncol = 3)
##
ggsave("./Fig5_TTrange_hist.pdf",hist_merge, width = 12.5, height= 12.5, units="in")
