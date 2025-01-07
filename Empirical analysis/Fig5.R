# setwd()
#---- Global option ----
clus_no = 4 # number of clusters
axis_size <- 14
text_size <- 12
title_size<- 16
col_list = c("#D05A6E", "#FB966E", "#51A8DD", "#5DAC81") # colors for cluster 1-4

#---- load data ----
library(dplyr)
df.all <- read.csv('./FAV_hv_data.csv',header = TRUE)
all_data_ = df.all %>% tidyr::drop_na(TTrange, W_length, Weight, B_length) # remove data with NA
all_data = all_data_ %>% dplyr::select(CTmin, CTmax, Family, Species, Location, Sample.Size, STmean, DTR, STR, Group) #select traits needed for analysis

#---- K-means cluster ----
library(factoextra)
k_clust = all_data_ %>% select(CTmin, CTmax)
set.seed(240802)
kmeans.cluster <- kmeans(k_clust, centers=clus_no) # clustering data into 4 clusters

### uniform cluster numbers (so that every result can be same color-code): 
for(i in c(1:length(kmeans.cluster$cluster))){
  if(kmeans.cluster$cluster[i] == 3) kmeans.cluster$cluster[i] = 1 # cluster 3 -> 1 (red)
  else if(kmeans.cluster$cluster[i] == 4) kmeans.cluster$cluster[i] = 3 # cluster 4 -> 3 (blue)
  else if(kmeans.cluster$cluster[i] == 1) kmeans.cluster$cluster[i] = 4 # cluster 1 -> 4 (green)
}
all_data$cluster = kmeans.cluster$cluster

## plot Fig.5A
p <- fviz_cluster(kmeans.cluster,
                  data = k_clust,
                  show.clust.cent = TRUE,
                  geom = c("point"),
                  frame.type = "norm")+
  coord_cartesian(xlim=c(-3,5), ylim = c(-6,3))+
  scale_colour_manual(values = scales::alpha(col_list,0.7))+
  scale_fill_manual(values = col_list)+ 
  scale_shape_manual(values = c(16,16,16,16))+
  theme_bw()+
  ggtitle('')+
  xlab('Critical thermal minimum (standardized)')+
  ylab('Critical thermal maximum (standardized)')+
  theme(axis.text = element_text(size = text_size),
        axis.title = element_text(size = axis_size),
        legend.position='inside',
        legend.position.inside=c(0.9,0.17),
        legend.background = element_blank())
p
# ggsave("./Fig5a.png",p, width = 5, height= 5, units="in")

#---- Cluster difference in groups ----
col_clus = paste0('clus',c(1:clus_no)) 
df_clus = as.data.frame(matrix(nrow = 9, ncol = clus_no+5))
colnames(df_clus) = c('Group', 'Location', 'Tmean', 'STR', 'DTR', col_clus)
df_clus$Group = c(c('M1', 'M2'), c('T1', 'T2', 'T3'), c('C1', 'C2', 'C3', 'C4'))
df_clus$Location = c(rep('Malaysia',2),rep('Taiwan',3),rep('China',4))
group_ls = c('M1', 'M2', 'T1', 'T2', 'T3', 'C1', 'C2', 'C3', 'C4')

for(i in c(1:9)){
  temp_vec = c()
  grp = group_ls[i]
  grp_df = all_data %>% filter(Group == grp)
  
  N = nrow(grp_df)
  ## Calculate cluster ratio
  clus_vec = c()
  for(k in c(1:clus_no)){
    clus_vec[k] = nrow(grp_df %>% filter(cluster == k))/N
  }
  df_clus[(df_clus$Group == grp),][(ncol(df_clus)-clus_no+1):ncol(df_clus)] = clus_vec
  
  ## Calculate average ambient temperature
  temp_vec[1] = mean(grp_df$STmean)
  temp_vec[2] = mean(grp_df$STR)
  temp_vec[3] = mean(grp_df$DTR)
  df_clus[(df_clus$Group == grp),][3:5] = temp_vec
}

#---- Plot: Fig.5B-D ----
#### Function: lmer regression model ####
library(lme4) #lmer
library(MuMIn) #r.squaredGLMM
library(ggeffects)
lmer_fun = function(df,temp_type,clus_num){
  coll_index = which(colnames(df) == 'Location')
  colx_index = which(colnames(df) == temp_type)
  coly_index = which(colnames(df) == clus_num)
  df_md = df[,c(coll_index,colx_index,coly_index)] # select specific temperature, cluster ratios, and location columns
  df_md = dplyr::rename(df_md, 'x' = temp_type, 'y' = clus_num)
  
  md_x <- lmer(y~ x+ (1|Location), data = df_md)
  result <- predict_response(md_x, 'x') # boundaries of confident interval
  sum=summary(md_x)
  coef_vec = c(round(sum$coefficients[2,1],2),round(r.squaredGLMM(md_x)[1],2),round(car::Anova(md_x)[['Pr(>Chisq)']],4)) # slope, r-square, p-value
  result_ls = list(result, coef_vec)
  
  # print(sum)
  # print(coef_vec)
  return(result_ls)
}

#### Function: ggplot ####
ggplot_fun = function(result_list, temp_type, clus_num, col, range.x, range.y){
  # adjust position of assemblage labels
  if(temp_type == 'Tmean') txt.mv.x = 1.7
  else if(temp_type == 'DTR') txt.mv.x = 0.5
  else if(temp_type == 'STR') txt.mv.x = 2
  
  df_pt = df_clus[,c(temp_type,clus_num,'Group')] # data-frame for plotting data-points
  colnames(df_pt) = c('x', 'y', 'Group')
  
  # coefficient text
  coef = result_list[[2]]
  slope = coef[1]
  r = coef[2]
  p = coef[3]
  
  if(p < 0.05 && p >= 0.01){
    coef_t = paste0('slope=',slope,', r2=',r,', p=',p,'*')
  }else if(p < 0.01 && p >= 0.001){
    coef_t =  paste0('slope=',slope,', r2=',r,', p=',p,'**')
  }else if(p < 0.001){
    coef_t =  paste0('slope=',slope,', r2=',r,', p<0.001***')
  }else{
    coef_t = paste0('slope=',slope,', r2=',r,', p=',p)
  }
  
  # plots
  df = result_list[[1]] # extract 95-c.i. df created in function lmer
  if(p > 0.05){ # if p> 0.05, dashed regression line
    assign('a', ggplot(df, aes(x = x, y = predicted))+
             geom_line(color=col, lty='dashed')+ # dashed regression line
             geom_point(aes(x = df_pt$x, y = df_pt$y), color=col)+
             geom_text(aes(label=df_pt$Group, x=df_pt$x+txt.mv.x, y=df_pt$y), color=col)+ # assemblage labels
             coord_cartesian(xlim=range.x, ylim=range.y)+
             xlab('')+
             ylab('')+
             theme_bw()+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none')+ 
             annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color='black')) # regression model coefficient
  }else{ # if p<= 0.05, solid regression line
    assign('a', ggplot(df, aes(x = x, y = predicted))+
             geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill=col, alpha = 0.1)+
             geom_line(color=col)+ # solid regression line
             geom_point(x = df_pt$x, y = df_pt$y, color=col)+
             geom_text(aes(label=df_pt$Group, x=df_pt$x+txt.mv.x, y=df_pt$y), color=col)+
             coord_cartesian(xlim=range.x, ylim=range.y)+
             xlab('')+
             ylab('')+
             theme_bw()+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none')+ 
             annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color='black'))
  }
  return(a)
}

#### Main plots Fig.5B-D ####
plot_ls_mean = list()
plot_ls_dtr = list()
plot_ls_str = list()
range.y=c(-0.01,1) # set y-axis limits
#--
for(temp_type in c('Tmean', 'DTR', 'STR')){
  for(i in c(1:clus_no)){
    if(temp_type == 'Tmean'){ ##-- 1. Tmean
      col=col_list[i]
      range.x=c(0,25) # set x-axis limits
      clus_num = paste0('clus',i)
      # function for coefficient
      result_list = lmer_fun(df_clus,temp_type,clus_num)
      # function for plotting
      assign(paste0(temp_type,i),ggplot_fun(result_list, temp_type, clus_num, col, range.x, range.y))
    }else if(temp_type == 'DTR'){ ##-- 2. DTR
      col=col_list[i]
      range.x=c(5,11)
      clus_num = paste0('clus',i)
      # function for coefficient
      result_list = lmer_fun(df_clus,temp_type,clus_num)
      # function for plotting
      assign(paste0(temp_type,i),ggplot_fun(result_list, temp_type, clus_num, col, range.x, range.y))
    }else if(temp_type == 'STR'){ ##-- 3. STR
      col=col_list[i]
      range.x=c(5,35)
      clus_num = paste0('clus',i)
      # function for coefficient
      result_list = lmer_fun(df_clus,temp_type,clus_num)
      # function for plotting
      assign(paste0(temp_type,i),ggplot_fun(result_list, temp_type, clus_num, col, range.x, range.y))
    }
  }
}
plot_ls_mean = mget(c(paste0('Tmean',c(1:clus_no))))
plot_ls_dtr = mget(c(paste0('DTR',c(1:clus_no))))
plot_ls_str = mget(c(paste0('STR',c(1:clus_no))))

#### Combined plots ####
library(ggpubr)
fig_mean<- ggarrange(plot_ls_mean[[1]],plot_ls_mean[[2]],plot_ls_mean[[3]],plot_ls_mean[[4]],nrow=1)
fig_comp1<- annotate_figure(fig_mean, bottom= text_grob("Average annual mean temperature (°C)", y = 1.5, size=axis_size))

fig_dtr<- ggarrange(plot_ls_dtr[[1]],plot_ls_dtr[[2]],plot_ls_dtr[[3]],plot_ls_dtr[[4]],nrow=1)
fig_comp2<- annotate_figure(fig_dtr, bottom= text_grob("Diurnal temperature range (°C)", y = 1.5, size=axis_size))

fig_str<- ggarrange(plot_ls_str[[1]],plot_ls_str[[2]],plot_ls_str[[3]],plot_ls_str[[4]],nrow=1)
fig_comp3<- annotate_figure(fig_str, bottom= text_grob("Seasonal temperature range (°C)", y = 1.5, size=axis_size))

fig_glmm<- ggarrange(fig_comp1,fig_comp2,fig_comp3,ncol=1,labels = c('B','C','D'))
fig_glmm<- annotate_figure(fig_glmm, left= text_grob("Cluster ratio", rot=90, size=axis_size))
ggsave(paste0("./Fig5B_D.pdf"),fig_glmm, width = 15, height= 12, units="in")

# fig_all<- ggarrange(p,fig_glmm,ncol=1,heights=c(1,3),labels = c('A',''))+ bgcolor("white")
# fig_all
# ggsave(paste0("./Fig5_full.pdf"),fig_all, width = 15, height= 20, units="in")
