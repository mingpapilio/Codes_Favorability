library(dplyr)
library(lme4)
library(MuMIn) 
library(AICcmodavg)

df.all <- read.csv('./FAV_hv_data.csv',header = TRUE)
all_data = df.all %>% dplyr::select(TTrange, Location, STmean, DTR, STR, Group) # select traits needed for analysis

#---- Global option ----
axis_size <- 14
text_size <- 12
title_size<- 16
theme_set(theme_bw())

#### Calculation Mean and SD of TTrange of each assemblage ####
df_TT = as.data.frame(matrix(nrow = 9, ncol = 8))
colnames(df_TT) = c('Location', 'Group', 'colors', 'TTMean', 'TTSD', 'Tmean', 'DTR', 'STR')
df_TT$Location = c(rep('Malaysia',2),rep('Taiwan',3),rep('China',4))
df_TT$Group = c(c('M1', 'M2'), c('T1', 'T2', 'T3'), c('C1', 'C2', 'C3', 'C4'))
df_TT$colors = c(rep('#ED784A',2),rep('#5DAC81',3),rep('#005CAF',4))
group_ls = c('M1', 'M2', 'T1', 'T2', 'T3', 'C1', 'C2', 'C3', 'C4')

## Calculate average ambient temperature of assemblages
Tmean_avg = all_data %>% group_by(Group) %>% summarise(Tmean = mean(STmean))
DTR_avg = all_data %>% group_by(Group) %>% summarise(DTR = mean(DTR))
STR_avg = all_data %>% group_by(Group) %>% summarise(STR = mean(STR))

## Calculate mean and SD of TTrange of assemblages and fill df_TT
for(i in c(1:9)){
  grp = group_ls[i]
  grp_df = all_data %>% filter(Group == grp)
  
  df_TT[(df_TT$Group == grp),]$TTMean = mean(grp_df$TTrange)
  df_TT[(df_TT$Group == grp),]$TTSD = sd(grp_df$TTrange)
  df_TT[(df_TT$Group == grp),]$Tmean = Tmean_avg[(Tmean_avg$Group == grp),]$Tmean
  df_TT[(df_TT$Group == grp),]$DTR = DTR_avg[(DTR_avg$Group == grp),]$DTR
  df_TT[(df_TT$Group == grp),]$STR = STR_avg[(STR_avg$Group == grp),]$STR
}

#### Fig.4D (Tmean) ####
####└　　GLMm model ####
md_x <- lmer(TTMean~ Tmean+ (1|Location), data = df_TT)
result_mean <- predict_response(md_x, "Tmean")
sum=summary(md_x)
coef_vec_mean = c(round(sum$coefficients[2,1],2),round(r.squaredGLMM(md_x)[1],2),round(car::Anova(md_x)[['Pr(>Chisq)']],4))

md_x <- lmer(TTSD~ Tmean+ (1|Location), data = df_TT)
result_sd <- predict_response(md_x, "Tmean")
sum=summary(md_x)
coef_vec_sd = c(round(sum$coefficients[2,1],2),round(r.squaredGLMM(md_x)[1],2),round(car::Anova(md_x)[['Pr(>Chisq)']],4))

####└　　ggplot ####
library(ggplot2)
library(ggeffects)
library(ggrepel)
col_list = c('#373C38','#373C38') # colors for 'Mean' and 'SD'
Mean_ls = list() # plot list for 'Tmean'

for(i in c(1:2)){
  nam <- paste("p", i, sep = "")
  if(i == 1){ # plot Mean of TTrange
    df = result_mean # GLMm model result
    df_pt = df_TT %>% select(Tmean, TTMean) %>% rename('y' = 'TTMean') # df of data-points
    coef = coef_vec_mean # coefficients of GLMm model
    label.y = 'Mean of TTrange' # label name of y-axis
    range.y = c(37,43) # y-axis limits
  }else if(i == 2){ # plot SD of TTrange
    df = result_sd
    df_pt = df_TT %>% select(Tmean, TTSD) %>% rename('y' = 'TTSD')
    coef = coef_vec_sd
    label.y = 'SD of TTrange'
    range.y = c(2,5)
  }
  ## coefficient text
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
  
  if(p > 0.05){ # no plotting ci area, and dashed regression line
    assign(nam, ggplot(df, aes(x = x, y = predicted, label = group_ls))+
             geom_line(color=col_list[i], lty='dashed')+
             geom_point(aes(x = df_pt$Tmean, y = df_pt$y), color=df_TT$colors)+
             coord_cartesian(xlim=c(3,21), ylim = range.y)+
             xlab('')+
             ylab(label.y)+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none',
                   panel.grid = element_line(linewidth = 0.1))+ 
             annotate(geom="text", x=15, y=range.y[2], label=coef_t, color=col_list[i])+
             geom_text_repel(x = df_pt$Tmean, y = df_pt$y, 
                             size = 5, color = col_list[i], vjust = 1.5))
  }else{ # plot ci area, and solid regression line
    assign(nam, ggplot(df, aes(x = x, y = predicted, label = group_ls))+
             geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill=col_list[i], alpha = 0.1)+
             geom_line(color=col_list[i])+
             geom_point(x = df_pt$Tmean, y = df_pt$y, color=df_TT$colors)+
             coord_cartesian(xlim=c(3,21), ylim = range.y)+
             xlab('')+
             ylab(label.y)+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none',
                   panel.grid = element_line(linewidth = 0.1))+ 
             annotate(geom="text", x=15, y=range.y[2], label=coef_t, color=col_list[i])+
             geom_text_repel(x = df_pt$Tmean, y = df_pt$y, 
                             size = 5, color = col_list[i], vjust = 1.5))
  }
}

Mean_ls <- mget(paste0("p", 1:2))

#### Fig.4G (DTR) ####
####└　　GLMm model ####
md_x <- lmer(TTMean~ DTR+ (1|Location), data = df_TT)
result_mean <- predict_response(md_x, "DTR")
sum=summary(md_x)
coef_vec_mean = c(round(sum$coefficients[2,1],2),round(r.squaredGLMM(md_x)[1],2),round(car::Anova(md_x)[['Pr(>Chisq)']],4))

md_x <- lmer(TTSD~ DTR+ (1|Location), data = df_TT)
result_sd <- predict_response(md_x, "DTR", interval='confidence', ci_level=0.95)
sum=summary(md_x)
coef_vec_sd = c(round(sum$coefficients[2,1],2),round(r.squaredGLMM(md_x)[1],2),round(car::Anova(md_x)[['Pr(>Chisq)']],4))

####└　　ggplot ####
col_list = c('#373C38','#373C38')
range.x = c(5,11)
dtr_ls = list()
for(i in c(1:2)){
  nam <- paste("p", i, sep = "")
  if(i == 1){
    df = result_mean
    df_pt = df_TT %>% select(DTR, TTMean) %>% rename('y' = 'TTMean')
    coef = coef_vec_mean
    label.y = 'Mean of TTrange'
    range.y = c(37,43)
  }else if(i == 2){
    df = result_sd
    df_pt = df_TT %>% select(DTR, TTSD) %>% rename('y' = 'TTSD')
    coef = coef_vec_sd
    label.y = 'SD of TTrange'
    range.y = c(2,5)
  }
  ## coefficient text
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
  
  if(p > 0.05){
    assign(nam, ggplot(df, aes(x = x, y = predicted, label = group_ls))+
             geom_line(color=col_list[i], lty='dashed')+
             geom_point(aes(x = df_pt$DTR, y = df_pt$y), color=df_TT$colors)+
             coord_cartesian(xlim=range.x, ylim = range.y)+
             xlab('')+
             ylab(label.y)+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none',
                   panel.grid = element_line(linewidth = 0.1))+ 
             annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color=col_list[i])+
             geom_text_repel(x = df_pt$DTR, y = df_pt$y, 
                             size = 5, color = col_list[i], vjust = 1.5))
  }else{
    assign(nam, ggplot(df, aes(x = x, y = predicted, label = group_ls))+
             geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill=col_list[i], alpha = 0.1)+
             geom_line(color=col_list[i])+
             geom_point(x = df_pt$DTR, y = df_pt$y, color=df_TT$colors)+
             coord_cartesian(xlim=range.x, ylim = range.y)+
             xlab('')+
             ylab(label.y)+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none',
                   panel.grid = element_line(linewidth = 0.1))+ 
             annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color=col_list[i])+
             geom_text_repel(x = df_pt$DTR, y = df_pt$y, 
                             size = 5, color = col_list[i], vjust = 1.5))
  }
}

dtr_ls <- mget(paste0("p", 1:2))

#### Fig.4J (STR) ####
####└　　GLMm model ####
md_x <- lmer(TTMean~ STR+ (1|Location), data = df_TT)
result_mean <- predict_response(md_x, "STR")
sum=summary(md_x)
coef_vec_mean = c(round(sum$coefficients[2,1],2),round(r.squaredGLMM(md_x)[1],2),round(car::Anova(md_x)[['Pr(>Chisq)']],4))

md_x <- lmer(TTSD~ STR+ (1|Location), data = df_TT)
result_sd <- predict_response(md_x, "STR")
sum=summary(md_x)
coef_vec_sd = c(round(sum$coefficients[2,1],2),round(r.squaredGLMM(md_x)[1],2),round(car::Anova(md_x)[['Pr(>Chisq)']],4))

####└　　ggplot ####
col_list = c('#373C38','#373C38')
range.x = c(9,32)
str_ls = list()
for(i in c(1:2)){
  nam <- paste("p", i, sep = "")
  if(i == 1){
    df = result_mean
    df_pt = df_TT %>% select(STR, TTMean) %>% rename('y' = 'TTMean')
    coef = coef_vec_mean
    label.y = 'Mean of TTrange'
    range.y = c(37,43)
  }else if(i == 2){
    df = result_sd
    df_pt = df_TT %>% select(STR, TTSD) %>% rename('y' = 'TTSD')
    coef = coef_vec_sd
    label.y = 'SD of TTrange'
    range.y = c(2,5)
  }
  ## coefficient text
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
  
  if(p > 0.05){
    assign(nam, ggplot(df, aes(x = x, y = predicted, label = group_ls))+
             geom_line(color=col_list[i], lty='dashed')+
             geom_point(aes(x = df_pt$STR, y = df_pt$y), color=df_TT$colors)+
             coord_cartesian(xlim=range.x, ylim = range.y)+
             xlab('')+
             ylab(label.y)+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none',
                   panel.grid = element_line(linewidth = 0.1))+ 
             annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color=col_list[i])+
             geom_text_repel(x = df_pt$STR, y = df_pt$y, 
                             size = 5, color = col_list[i], vjust = 1.5))
  }else{
    assign(nam, ggplot(df, aes(x = x, y = predicted, label = group_ls))+
             geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill=col_list[i], alpha = 0.1)+
             geom_line(color=col_list[i])+
             geom_point(x = df_pt$STR, y = df_pt$y, color=df_TT$colors)+
             coord_cartesian(xlim=range.x, ylim = range.y)+
             xlab('')+
             ylab(label.y)+
             theme(axis.text = element_text(size = text_size),
                   axis.title = element_text(size = axis_size),
                   legend.position = 'none',
                   panel.grid = element_line(linewidth = 0.1))+ 
             annotate(geom="text", x=mean(range.x), y=range.y[2], label=coef_t, color=col_list[i])+
             geom_text_repel(x = df_pt$STR, y = df_pt$y, 
                             size = 5, color = col_list[i], vjust = 1.5))
  }
}

str_ls <- mget(paste0("p", 1:2))

#### Combined plots ####
library(ggpubr)
fig_Tmean<- ggarrange(Mean_ls[[1]],Mean_ls[[2]],ncol=1,labels = c("D",""))
fig_comp_mean<- annotate_figure(fig_Tmean, 
                                bottom= text_grob("Average annual mean temperature (°C)", 
                                                  y = 1.5, size=axis_size))
fig_DTR<- ggarrange(dtr_ls[[1]],dtr_ls[[2]],ncol=1,labels = c("G",""))
fig_comp_dtr<- annotate_figure(fig_DTR, 
                               bottom= text_grob("Diurnal temperature range (°C)", 
                                                 y = 1.5, size=axis_size))
fig_STR<- ggarrange(str_ls[[1]],str_ls[[2]],ncol=1,labels = c("J",""))
fig_comp_str<- annotate_figure(fig_STR, 
                               bottom= text_grob("Seasonal temperature range (°C)", 
                                                 y = 1.5, size=axis_size))
blank = ggplot() + theme_void() # Add blank plot so the plots can be composed easier
fig_merge<- ggarrange(fig_comp_mean,fig_comp_dtr,fig_comp_str,blank,ncol=2,nrow=2)+ bgcolor("white")
fig_merge
ggsave("./Fig4DGJ.pdf",fig_merge, width = 8, height= 10, units="in")
