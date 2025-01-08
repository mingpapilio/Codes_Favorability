# setwd()
#----- global option ------
random_time = 100

########## Function ##########
#---- p-value function (lm) ----
p_value_lm = function(df_x, df_y, xlim, ylim){

  md <- lm(df_y~ df_x, data = df_)
  sum = summary(md)
  slope = round(sum[["coefficients"]][2], 2)
  r = round(r.squaredGLMM(md)[1], 2)
  p = round(sum[["coefficients"]][8], 3)

  text_pos_x = sum(xlim)/2

  if(p < 0.05 && p >= 0.01){
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,'*',sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,'**',sep=''), cex=1)
  }else if(p < 0.001){
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p < 0.001***',sep=''), cex=1)
  }else{
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,',p = ',p,sep=''), cex=1)
  }
  return(p)
}

########## Main ##########
#---- load data ----
df.all <- read.csv('./FAV_hv_data.csv',header = TRUE)
all_data = df.all %>% dplyr::select(CTmin, CTmax, Family, Species, Location, Sample.Size, STmean, DTR, STR, Group) # select traits needed for analysis

# Standardize traits
library(vegan)
all_data$CTmax_sd = decostand(all_data$CTmax,"standardize",na.rm=T)[,1]
all_data$CTmin_sd = decostand(all_data$CTmin,"standardize",na.rm=T)[,1]

#---- Calculate hyper-volume (each assemblage) ----
# create result data frame
df_ = as.data.frame(matrix(nrow = 9, ncol = 3))
colnames(df_) = c('Location', 'group', 'Thermal_hv')
df_$Location = c(rep('Malaysia',2),rep('Taiwan',3),rep('China',4))
df_$group = c(c('M1', 'M2'), c('T1', 'T2', 'T3'), c('C1', 'C2', 'C3', 'C4'))

# hyper-volume bootstrap
library(hypervolume)
for(loc in c('Malaysia', 'Taiwan', 'China')){
  if(loc == 'Malaysia') group = c('M1', 'M2')
  else if(loc == 'Taiwan') group = c('T1', 'T2', 'T3')
  else if(loc == 'China') group = c('C1', 'C2', 'C3', 'C4')
  
  for(num in group){
    df <- all_data %>% filter(Location == loc & Group == num)
    df_random_thermal_hv = c()
    
    i = 1
    
    for(i in seq(1,random_time)){
      set.seed(20231101+i)
      df2 = sample_n(df, 20) # randomly pick 20 species in each group to calculate hyper-volume
      
      # create weighting data-frame (sample size of each species) for each group 
      comm1 = data.frame(row.names = df2$Species)
      comm1$n = NA
      comm1$n[match(df2$Species, rownames(comm1))] <- df2$Sample.Size
      comm1[is.na(comm1)] <- 0
      comm1 = t(comm1)
      weight = comm1/sum(comm1)
      
      # calculating hyper-volume
      trait_a = dplyr::select(df2, CTmax_sd, CTmin_sd)
      rownames(trait_a) = df2$Species
      hvlist_a <- hypervolume::hypervolume_gaussian(trait_a, weight = weight)
      
      df_random_thermal_hv[i] = get_volume(hvlist_a)
      
      print(paste(i,"/",random_time,'/',loc,'/',num))
    }
    
    df_[(df_$Location == loc & df_$group == num),]$Thermal_hv = mean(df_random_thermal_hv)
  }
}

#---- Calculate mean N/Tmean/DTR/STR of each assemblage and add into hyper-volume dataframe ----
# Average N of species/Tmean/DTR/STR of each assemblage
N_avg = all_data %>% group_by(Group) %>% summarise(N = n())
Tmean_avg = all_data %>% group_by(Group) %>% summarise(Tmean = mean(STmean))
DTR_avg = all_data %>% group_by(Group) %>% summarise(DTR = mean(DTR))
STR_avg = all_data %>% group_by(Group) %>% summarise(STR = mean(STR))

# Adding temperature data in hyper-volume data-frame
df_$N_all = NA
df_$Tmean_all = NA
df_$DTR_all = NA
df_$STR_all = NA
for(loc in c('Malaysia','Taiwan','China')){
  if(loc == 'Malaysia') group = c('M1', 'M2')
  else if(loc == 'Taiwan') group = c('T1', 'T2', 'T3')
  else if(loc == 'China') group = c('C1', 'C2', 'C3', 'C4')
  
  for(num in group){
    df_[(df_$Location == loc & df_$group == num),]$N_all = N_avg[(N_avg$Location == loc & Tmean_avg$Group == num),]$N
    df_[(df_$Location == loc & df_$group == num),]$Tmean_all = Tmean_avg[(Tmean_avg$Location == loc & Tmean_avg$Group == num),]$Tmean
    df_[(df_$Location == loc & df_$group == num),]$DTR_all = DTR_avg[(DTR_avg$Location == loc & DTR_avg$Group == num),]$DTR
    df_[(df_$Location == loc & df_$group == num),]$STR_all = STR_avg[(STR_avg$Location == loc & STR_avg$Group == num),]$STR
  }
}

# write.csv(df_, file = './hv_result.csv', row.names = F)

#---- Plot GLMm (Fig.4B,4E,4H) ----
# df_ = read.csv('./hv_result.csv')

pdf(file = './Fig4BEH.pdf', width = 15, height = 5, onefile = TRUE)

par(mfrow=c(3,1))

## Fig.4B (Tmean)
df_x = df_$Tmean_all
df_y = df_$Thermal_hv

plot(x=df_x, y=df_y, las=1, cex= 2, 
     xlab = '', ylab = '', main = 'STmean',
     xlim = c(0, 25), ylim = c(5,25), type="n", axes = F)
box()

axis(side=1, at=seq(0, 25, 5),cex.axis=1.5)
mtext(side=1, line=3, 'Averaged annual mean temperature (°C)', col="black", font=1, cex=1)
axis(side=2, at=seq(5,25,5),cex.axis=1.5, las=1)
mtext(side=2, line=3, 'Hypervolume of thermal traits', col="black", font=1, cex=1)

# GLMm model
md_x <- lm(df_y~ df_x, data = df_)
mydata.x <- data.frame(expand.grid(df_x = seq(0, 25, 0.1)))
conf_interval <- as.data.frame(predict(md_x, newdata=data.frame(x=mydata.x), interval="confidence",
                                       level = 0.95))
p_ = p_value_lm(df_x, df_y, c(0, 25), 5) # function for producing coefficients
colors <- c("#005CAF", "#ED784A", "#5DAC81")
points(x=df_x, y=df_y,
       col=colors[factor(df_$Location)], pch=16, cex=1.2) # points are colored by locations
text(df_y~ df_x, labels=df_$group, data=df_, cex=1, font=1, pos=4) # assemblage IDs besides points
clip(min(df_x)-0.5,max(df_x)+0.5, 4.5, 100) # set limitation on regression line and ci area
if(p_ > 0.05) lty = 2 else{ # draw ci interval if p<=0.05, else skip
  lty = 1
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(conf_interval$lwr,rev(conf_interval$upr)),
          col=scales::alpha('#787878', 0.2),border=NA)
}
abline(md_x,col='black',lwd=1.5,lty=lty) # draw regression line

## Fig.4E (DTR)
df_x = df_$DTR_all
df_y = df_$Thermal_hv

plot(x=df_x, y=df_y, las=1, cex= 2, 
     xlab = '', ylab = '', main = 'DTR',
     xlim = c(5, 11), ylim = c(5,25), type="n", axes = F)
box()

axis(side=1, at=seq(5, 11, 1),cex.axis=1.5)
mtext(side=1, line=3, 'Diurnal temperature range (°C)', col="black", font=1, cex=1)
axis(side=2, at=seq(5,25, 5),cex.axis=1.5, las=1)
mtext(side=2, line=3, 'Hypervolume of thermal traits', col="black", font=1, cex=1)

md_x <- lm(df_y~ df_x, data = df_)
mydata.x <- data.frame(expand.grid(df_x = seq(5, 10, 0.1)))
conf_interval <- as.data.frame(predict(md_x, newdata=data.frame(x=mydata.x), interval="confidence",
                                       level = 0.95))
p_ = p_value_lm(df_x, df_y, c(5, 10), 5)
colors <- c("#005CAF", "#ED784A", "#5DAC81")
points(x=df_x, y=df_y,
       col=colors[factor(df_$Location)], pch=16, cex=1.2)
text(df_y~ df_x, labels=df_$group, data=df_, cex=1, font=1, pos=4)
clip(min(df_x)-0.5,max(df_x)+0.5, 4.5, 100)
if(p_ > 0.05) lty = 2 else{ 
  lty = 1
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(conf_interval$lwr,rev(conf_interval$upr)),
          col=scales::alpha('#787878', 0.2),border=NA)
}
abline(md_x,col='black',lwd=1.5,lty=lty)

## Fig.4H (STR)
df_x = df_$STR_all
df_y = df_$Thermal_hv

plot(x=df_x, y=df_y, las=1, cex= 2, 
     xlab = '', ylab = '', main = 'STR',
     xlim = c(5, 40), ylim = c(5,25), type="n", axes = F)
box()

axis(side=1, at=seq(5, 40, 5),cex.axis=1.5)
mtext(side=1, line=3, 'Seasonal temperature range (°C)', col="black", font=1, cex=1)
axis(side=2, at=seq(5,25,5),cex.axis=1.5, las=1)
mtext(side=2, line=3, 'Hypervolume of thermal traits', col="black", font=1, cex=1)

md_x <- lm(df_y~ df_x, data = df_)
mydata.x <- data.frame(expand.grid(df_x = seq(5, 40, 0.1)))
conf_interval <- as.data.frame(predict(md_x, newdata=data.frame(x=mydata.x), interval="confidence",
                                       level = 0.95))
p_ = p_value_lm(df_x, df_y, c(5, 40), 5)
colors <- c("#005CAF", "#ED784A", "#5DAC81")
points(x=df_x, y=df_y,
       col=colors[factor(df_$Location)], pch=16, cex=1.2)
text(df_y~ df_x, labels=df_$group, data=df_, cex=1, font=1, pos=4)
clip(min(df_x)-0.5,max(df_x)+0.5, 4.5, 100)
if(p_ > 0.05) lty = 2 else{ 
  lty = 1
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(conf_interval$lwr,rev(conf_interval$upr)),
          col=scales::alpha('#787878', 0.2),border=NA)
}
abline(md_x,col='black',lwd=1.5,lty=lty)

dev.off()