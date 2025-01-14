############# Use same data in Fig.2 #############
#--
library(dplyr)
library(vegan)
library(MuMIn) #r.squaredGLMM
library(nlme) #lme
library(AICcmodavg)
library(data.table) #fread
library(factoextra) #fviz_cluster

#---- Global options ----
clus_no = 4
axis_size<- 14
text_size<- 12
title_size<- 16
theme_set(theme_bw())
if(clus_no != 2) col_list = c("#D05A6E", "#FB966E", "#51A8DD", "#5DAC81")else col_list = c("#D05A6E", "#51A8DD")

#---- Fig. K-Mean cluster + GLMm ----
#### └ load data ####
spcs_data_mean <- fread("./mean_svar2.5_lvar0.0/species.txt")
spcs_data_svar <- fread("./svar_mean18.0_lvar0.0/species.txt")
spcs_data_lvar <- fread("./lvar_mean18.0_svar2.5/species.txt")

#### └ K-means cluster ####
#### └　 Fig.3a. Mean ambient temperature ####
k_clust = spcs_data_mean %>% select(CTmin, CTmax)

set.seed(2408021)
kmeans.cluster <- kmeans(k_clust, centers=clus_no)
### uniform cluster numbers: 
## 4 clusters
for(i in c(1:length(kmeans.cluster$cluster))){
  if(kmeans.cluster$cluster[i] == 4) kmeans.cluster$cluster[i] = 1
  else if(kmeans.cluster$cluster[i] == 1) kmeans.cluster$cluster[i] = 2
  else if(kmeans.cluster$cluster[i] == 2) kmeans.cluster$cluster[i] = 3
  else if(kmeans.cluster$cluster[i] == 3) kmeans.cluster$cluster[i] = 4
}
spcs_data_mean$cluster = kmeans.cluster$cluster

## plot
a1 <- fviz_cluster(kmeans.cluster,
                  data = k_clust,
                  show.clust.cent = TRUE,
                  geom = c("point"),
                  frame.type = "norm")+
  coord_cartesian(xlim=c(-6,4), ylim = c(-7,7))+
  scale_colour_manual(values = scales::alpha(col_list,0.3))+
  scale_fill_manual(values = col_list)+ 
  scale_shape_manual(values = c(16,16,16,16))+
  theme_bw()+
  ggtitle('')+
  xlab('Critical thermal minimum (°C)')+
  ylab('Critical thermal maximum (°C)')+
  theme(axis.text = element_text(size = text_size),
        axis.title = element_text(size = axis_size),
        legend.position=c(0.9,0.17), legend.background = element_blank())
a1

#### └　 Fig.3c. Short-term ####
k_clust = spcs_data_svar %>% select(CTmin, CTmax)

set.seed(2408022)
kmeans.cluster <- kmeans(k_clust, centers=clus_no)
### uniform cluster numbers: 
## 4 clusters
for(i in c(1:length(kmeans.cluster$cluster))){
  if(kmeans.cluster$cluster[i] == 2) kmeans.cluster$cluster[i] = 1
  else if(kmeans.cluster$cluster[i] == 4) kmeans.cluster$cluster[i] = 2
  else if(kmeans.cluster$cluster[i] == 1) kmeans.cluster$cluster[i] = 4
}
spcs_data_svar$cluster = kmeans.cluster$cluster

## plot
b1 <- fviz_cluster(kmeans.cluster,
                  data = k_clust,
                  show.clust.cent = TRUE,
                  geom = c("point"),
                  frame.type = "norm")+
  coord_cartesian(xlim=c(-6,4), ylim = c(-7,7))+
  scale_colour_manual(values = scales::alpha(col_list,0.3))+
  scale_fill_manual(values = col_list)+ 
  scale_shape_manual(values = c(16,16,16,16))+
  theme_bw()+
  ggtitle('')+
  xlab('Critical thermal minimum (°C)')+
  ylab('Critical thermal maximum (°C)')+
  theme(axis.text = element_text(size = text_size),
        axis.title = element_text(size = axis_size),
        legend.position=c(0.9,0.17), legend.background = element_blank())
  #legend.position='none'
b1

#### └　 Fig.3e. Long-term ####
k_clust = spcs_data_lvar %>% select(CTmin, CTmax)

set.seed(2408023) # seeds: 4clus:2408023, 3clus:332408023, 2clus:222408023
kmeans.cluster <- kmeans(k_clust, centers=clus_no)
### uniform cluster numbers: 
## 4 clusters
for(i in c(1:length(kmeans.cluster$cluster))){
  if(kmeans.cluster$cluster[i] == 4) kmeans.cluster$cluster[i] = 1
  else if(kmeans.cluster$cluster[i] == 3) kmeans.cluster$cluster[i] = 2
  else if(kmeans.cluster$cluster[i] == 2) kmeans.cluster$cluster[i] = 3
  else if(kmeans.cluster$cluster[i] == 1) kmeans.cluster$cluster[i] = 4
}
spcs_data_lvar$cluster = kmeans.cluster$cluster

## plot  
c1 <- fviz_cluster(kmeans.cluster,
                  data = k_clust,
                  show.clust.cent = TRUE,
                  geom = c("point"),
                  frame.type = "norm")+
  coord_cartesian(xlim=c(-6,4), ylim = c(-7,7))+
  scale_colour_manual(values = scales::alpha(col_list,0.3))+
  scale_fill_manual(values = col_list)+ 
  scale_shape_manual(values = c(16,16,16,16))+
  theme_bw()+
  ggtitle('')+
  xlab('Critical thermal minimum (°C)')+
  ylab('Critical thermal maximum (°C)')+
  theme(axis.text = element_text(size = text_size),
        axis.title = element_text(size = axis_size),
        legend.position=c(0.9,0.17), legend.background = element_blank())
c1

#### └ Cluster ratio in groups plot (Fig.3b,3d,3f) ####
range.y = c(0,0.5)
label.x = ''
label.y = ''

#### └　 Function ####
fun_glmm_plot = function(df_clus, temp_type, range.x, range.y, label.x, label.y){
  for(i in c(1:clus_no)){
    df = df_clus    
    fea_type = paste0('clus',i)
    nam = paste0("p", i)
    x_pos = which(colnames(df) == temp_type)
    y_pos = which(colnames(df) == fea_type)
    colnames(df)[x_pos] = 'Group'
    colnames(df)[y_pos] = 'y'
    col = col_list[i]
    
    ## coefficient text
    md <- lm(y~ Group, data = df)
    sum = summary(md)
    slope = round(sum[["coefficients"]][2], 2)
    r = round(r.squaredGLMM(md)[1], 2)
    p = round(sum[["coefficients"]][8], 3)
    
    if(p < 0.05 && p >= 0.01){
      coef_t = paste0('slope=',slope,', r2=',r,', p=',p,'*')
    }else if(p < 0.01 && p >= 0.001){
      coef_t =  paste0('slope=',slope,', r2=',r,', p=',p,'**')
    }else if(p < 0.001){
      coef_t =  paste0('slope=',slope,', r2=',r,', p<0.001***')
    }else{
      coef_t = paste0('slope=',slope,', r2=',r,', p=',p)
    }
    
    if(p > 0.05){
      assign(nam, ggplot(df, aes(x = Group, y = y))+
               geom_smooth(method = "lm", se = FALSE, lty='dashed', 
                           color=col, alpha = 0.3) + 
               geom_point(color=col,alpha=0.5,size=2.5)+
               theme(axis.text = element_text(size = text_size),
                     legend.position = 'none')+
               coord_cartesian(xlim=range.x, ylim = range.y)+
               xlab(label.x)+
               ylab(label.y)+
               theme(axis.text = element_text(size = text_size),
                     axis.title = element_text(size = axis_size),
                     legend.position = 'none')+ 
               theme_bw()+
               annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color=col))
    }else{
      assign(nam, ggplot(df, aes(x = Group, y = y))+
               geom_smooth(method = "lm", se = TRUE, color=col, fill = col, alpha = 0.3) + 
               geom_point(color=col,alpha=0.5,size=2.5)+
               theme(axis.text = element_text(size = text_size),
                     legend.position = 'none')+
               coord_cartesian(xlim=range.x, ylim = range.y)+
               xlab(label.x)+
               ylab(label.y)+
               theme(axis.text = element_text(size = text_size),
                     axis.title = element_text(size = axis_size),
                     legend.position = 'none')+ 
               theme_bw()+
               annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color=col))
    }
  }
  plot_list <- mget(paste0("p", 1:clus_no))
  return(plot_list)
}

## lm tables ##
fun_lm_table = function(df_clus, temp_type){
  for(i in c(1:clus_no)){
    df = df_clus    
    fea_type = paste0('clus',i)
    nam = paste0("p", i)
    x_pos = which(colnames(df) == temp_type)
    y_pos = which(colnames(df) == fea_type)
    colnames(df)[x_pos] = 'Group'
    colnames(df)[y_pos] = 'y'
    col = col_list[i]
    
    ## coefficient text
    md <- lm(y~ Group, data = df)
    print(paste('Cluster', i))
    print(summary(md))
    print(round(r.squaredGLMM(md)[1], 4))
  } 
}

#### └　 Fig.3b. Mean ambient temperature ####
col_clus = paste0('clus',c(1:clus_no)) 
df_clus_mean = as.data.frame(matrix(nrow = 11, ncol = clus_no+1))
colnames(df_clus_mean) = c('Tmean', col_clus)
df_clus_mean$Tmean = unique(spcs_data_mean$Tmean)

for(i in c(1:11)){
  clus_vec = c()
  grp = unique(spcs_data_mean$Tmean)[i]
  grp_df = spcs_data_mean %>% filter(Tmean == grp)
  
  ## Calculate cluster ratio
  N = nrow(grp_df)
  for(k in c(1:clus_no)){
    clus_vec[k] = nrow(grp_df %>% filter(cluster == k))/N
  }
  df_clus_mean[(df_clus_mean$Tmean == grp),][(ncol(df_clus_mean)-clus_no+1):ncol(df_clus_mean)] = clus_vec
}

## plot glmm ##
range.x=c(14,22)
plot_ls_mean = fun_glmm_plot(df_clus_mean, 'Tmean', range.x, range.y, label.x, label.y)
## Print GLMm coefficient table ##
fun_lm_table(df_clus_mean, 'Tmean')

#### └　 Fig.3d. Short-term ####
col_clus = paste0('clus',c(1:clus_no)) 
df_clus_svar = as.data.frame(matrix(nrow = 11, ncol = clus_no+1))
colnames(df_clus_svar) = c('Tsd_short', col_clus)
df_clus_svar$Tsd_short = unique(spcs_data_svar$Tsd_short)

for(i in c(1:11)){
  clus_vec = c()
  grp = unique(spcs_data_svar$Tsd_short)[i]
  grp_df = spcs_data_svar %>% filter(Tsd_short == grp)
  
  ## Calculate cluster ratio
  N = nrow(grp_df)
  for(k in c(1:clus_no)){
    clus_vec[k] = nrow(grp_df %>% filter(cluster == k))/N
  }
  df_clus_svar[(df_clus_svar$Tsd_short == grp),][(ncol(df_clus_svar)-clus_no+1):ncol(df_clus_svar)] = clus_vec
}

## plot glmm ##
range.x=c(0,5)
plot_ls_svar = fun_glmm_plot(df_clus_svar, 'Tsd_short', range.x, range.y, label.x, label.y)
## Print GLMm coefficient table ##
fun_lm_table(df_clus_svar, 'Tsd_short')

#### └　 Fig.3f. Long-term ####
col_clus = paste0('clus',c(1:clus_no)) 
df_clus_lvar = as.data.frame(matrix(nrow = 11, ncol = clus_no+1))
colnames(df_clus_lvar) = c('Tsd', col_clus)
df_clus_lvar$Tsd = unique(spcs_data_lvar$Tsd)

for(i in c(1:11)){
  clus_vec = c()
  grp = unique(spcs_data_lvar$Tsd)[i]
  grp_df = spcs_data_lvar %>% filter(Tsd == grp)
  
  ## Calculate cluster ratio
  N = nrow(grp_df)
  for(k in c(1:clus_no)){
    clus_vec[k] = nrow(grp_df %>% filter(cluster == k))/N
  }
  df_clus_lvar[(df_clus_lvar$Tsd == grp),][(ncol(df_clus_lvar)-clus_no+1):ncol(df_clus_lvar)] = clus_vec
}

## plot glmm ##
range.x=c(0,5)
plot_ls_lvar = fun_glmm_plot(df_clus_lvar, 'Tsd', range.x, range.y, label.x, label.y)
## Print GLMm coefficient table ##
fun_lm_table(df_clus_lvar, 'Tsd')

#### └ Combined plots ####
library(ggpubr)
widths_ = c(1,1)
plot_wid = 14
plot_hgt = 18
col_ = 2
row_ = 2

a3<- ggarrange(plot_ls_mean[[1]],plot_ls_mean[[2]],plot_ls_mean[[3]],plot_ls_mean[[4]],ncol=col_,nrow=row_) #,plot_ls_mean[[3]],plot_ls_mean[[4]]
a3 <- annotate_figure(a3, bottom= text_grob("Mean ambient temperature (°C)", size=axis_size, y=1))
fig_comp1 <- ggarrange(a1,a3,nrow=1,widths=widths_,labels = c("A","B"))
fig_comp1 <- annotate_figure(fig_comp1, top= text_grob("Mean ambient temperature", size=title_size, face = "bold"))

b3<- ggarrange(plot_ls_svar[[1]],plot_ls_svar[[2]],plot_ls_svar[[3]],plot_ls_svar[[4]],ncol=col_,nrow=row_) #,plot_ls_svar[[3]],plot_ls_svar[[4]]
b3 <- annotate_figure(b3, bottom= text_grob("Short-term variability (°C)", size=axis_size, y=1))
fig_comp2 <- ggarrange(b1,b3,nrow=1,widths=widths_,labels = c("C","D"))
fig_comp2 <- annotate_figure(fig_comp2, top= text_grob("Short-term variability", size=title_size, face = "bold"))

c3<- ggarrange(plot_ls_lvar[[1]],plot_ls_lvar[[2]],plot_ls_lvar[[3]],plot_ls_lvar[[4]],ncol=col_,nrow=row_) #,plot_ls_lvar[[3]],plot_ls_lvar[[4]]
c3 <- annotate_figure(c3, bottom= text_grob("Long-term variability (°C)", size=axis_size, y=1))
fig_comp3 <- ggarrange(c1,c3,nrow=1,widths=widths_,labels = c("E","F"))
fig_comp3 <- annotate_figure(fig_comp3, top= text_grob("Long-term variability", size=title_size, face = "bold"))

fig_merge<- ggarrange(fig_comp1,fig_comp2,fig_comp3, ncol=1)+ bgcolor("white")
fig_merge
ggsave('./Fig3.pdf',fig_merge, width = plot_wid, height= plot_hgt, units="in")
