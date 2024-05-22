library(dplyr)
library(hypervolume)
library(vegan)
library(lme4) #lmer
library(MuMIn) #r.squaredGLMM
library(nlme) #lme
library(AICcmodavg)

#--------------
random_time = 100

########## Function ##########
#---- Plotting Function ----
plot_fun <- function(df_x, df_y, xlab, ylab, adj_x = 3, adj_y = 3){
  
  plot(x=df_x, y=df_y, las=1,cex= 1.5, 
       xlab = xlab, ylab = ylab,
       xlim=c(min(df_x), max(df_x)), ylim=c(min(df_y), max(df_y)), type="n")
  
  md_x <- lmer(df_y~ df_x+ (1|df_[['Location']]))
  mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x)-1, max(df_x)+1, 0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  
  md <- lme(df_y~ df_x, random= ~1|Location, data = df_, na.action=na.exclude)
  sum = summary(md)
  r = round(r.squaredGLMM(md)[1], 4)
  p = round(sum[["tTable"]][10], 4)
  
  if(p < 0.05 && p >= 0.01){
    lines(x=mydata.x$df_x,y=pred.x$fit,col="black",lwd=1,lty=1)
    text(max(df_x)-adj_x, min(df_y)-adj_y, paste('p = ',p,'*, r-square = ',r,sep=''), cex=0.8)
  }else if(p < 0.01 && p >= 0.001){
    lines(x=mydata.x$df_x,y=pred.x$fit,col="black",lwd=1,lty=1)
    text(max(df_x)-adj_x, min(df_y)-adj_y, paste('p = ',p,'**, r-square = ',r,sep=''), cex=0.8)
  }else if(p < 0.001){
    lines(x=mydata.x$df_x,y=pred.x$fit,col="black",lwd=1,lty=1)
    text(max(df_x)-adj_x, min(df_y)-adj_y, 'p < 0.001***, r-square = ',r, cex=0.8)
  }else{
    lines(x=mydata.x$df_x,y=pred.x$fit,col="black",lwd=1,lty=2)
    text(max(df_x)-adj_x, min(df_y)-adj_y, paste('p = ',p,sep=''), cex=0.8)
  }
  
  colors <- c("#005CAF","#ED784A","#5DAC81")
  points(x=df_x, y=df_y,
         col=colors[factor(df_$Location)], pch=16, cex=1)
  text(df_y~ df_x, labels=df_$group, data=df_, cex=0.6, font=2, pos=4)
  
}

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
all_data_ = df.all %>% tidyr::drop_na(CTrange, W_length, Weight, B_length)
all_data = all_data_ %>% dplyr::select(CTmin, CTmax, Family, Species, Location, Sample.Size, STmean, DTR, STR, Group) #Tmean

## grouping species into high, mid and low STmean, DTR and STR ####
all_data = all_data %>% mutate(Tmean_g = ifelse(Group == 'M1'|Group == 'M2'|Group == 'T1','High', 
                                                ifelse(Group == 'T2'|Group == 'C1'|Group == 'T3', 'Mid', 'Low')))
all_data = all_data %>% mutate(DTR_g = ifelse(Group == 'C1'|Group == 'C3'|Group == 'T1', 'High',
                                                  ifelse(Group == 'C2'|Group == 'C4'|Group == 'T2', 'Mid','Low')))
all_data = all_data %>% mutate(STR_g = ifelse(Group == 'C4'|Group == 'C3'|Group == 'C2', 'High',
                                                  ifelse(Group == 'C1'|Group == 'T2'|Group == 'T3', 'Mid','Low')))

# Standardize traits
all_data$CTmax_sd = decostand(all_data$CTmax,"standardize",na.rm=T)[,1]
all_data$CTmin_sd = decostand(all_data$CTmin,"standardize",na.rm=T)[,1]

#---- Calculate total hyper-volume (each assemblage) ----
# create result data frame
df_ = as.data.frame(matrix(nrow = 9, ncol = 3))
colnames(df_) = c('Location', 'group', 'Thermal_hv')
df_$Location = c(rep('Malaysia',2),rep('Taiwan',3),rep('China',4))
df_$group = c(c('M1', 'M2'), c('T1', 'T2', 'T3'), c('C1', 'C2', 'C3', 'C4'))

# hyper-volume bootstrap
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
# Average Tmean/DTR of 3 locations
N_avg = all_data %>% group_by(Location, Group) %>% summarise(N = n())
Tmean_avg = all_data %>% group_by(Location, Group) %>% summarise(Tmean = mean(STmean)) #%>% arrange(desc(Tmean))
DTR_avg = all_data %>% group_by(Location, Group) %>% summarise(DTR = mean(DTR)) #%>% arrange(desc(DTR))
STR_avg = all_data %>% group_by(Location, Group) %>% summarise(STR = mean(STR)) #%>% arrange(desc(STR))

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

df_$group_num = df_$group
df_$group = c('M1', 'M2', 'T1', 'T2', 'T3', 'C1', 'C2', 'C3', 'C4')

#---- Plot lm (Fig 4B-E) ----
pdf(file = './Fig4B-E.pdf', width = 10, height = 10, onefile = TRUE)

par(mfrow=c(2,2))

## N of species
df_x = df_$N_all
df_y = df_$Thermal_hv

plot(x=df_x, y=df_y, las=1, cex= 2, 
     xlab = '', ylab = '', main = 'N of species',
     xlim = c(20, 140), ylim = c(5,25), type="n", axes = F)
box()

axis(side=1, at=seq(20, 140, 20),cex.axis=1.5)
mtext(side=1, line=3, 'Number of species', col="black", font=1, cex=1)
axis(side=2, at=seq(5,25,5),cex.axis=1.5, las=1)
mtext(side=2, line=3, 'Hypervolume of thermal traits', col="black", font=1, cex=1)

md_x <- lm(df_y~ df_x, data = df_)
mydata.x <- data.frame(expand.grid(df_x = seq(20, 140, 0.1)))
conf_interval <- as.data.frame(predict(md_x, newdata=data.frame(x=mydata.x), interval="confidence",
                                       level = 0.95))
p_ = p_value_lm(df_x, df_y, c(20, 140), 5)
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

## Tmean
# df_x = df_$Tmean
df_x = df_$Tmean_all
df_y = df_$Thermal_hv

plot(x=df_x, y=df_y, las=1, cex= 2, 
     xlab = '', ylab = '', main = 'STmean',
     xlim = c(0, 25), ylim = c(5,25), type="n", axes = F)
box()

axis(side=1, at=seq(0, 25, 5),cex.axis=1.5)
mtext(side=1, line=3, 'Average annual mean temperature (°C)', col="black", font=1, cex=1)
axis(side=2, at=seq(5,25,5),cex.axis=1.5, las=1)
mtext(side=2, line=3, 'Hypervolume of thermal traits', col="black", font=1, cex=1)

md_x <- lm(df_y~ df_x, data = df_)
mydata.x <- data.frame(expand.grid(df_x = seq(0, 25, 0.1)))
conf_interval <- as.data.frame(predict(md_x, newdata=data.frame(x=mydata.x), interval="confidence",
                                       level = 0.95))
p_ = p_value_lm(df_x, df_y, c(0, 25), 5)
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

# DTR
# df_x = df_$DTR
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

# STR
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


#---- Fig 4F-H. 0.95 kernel density (high-mid-low) ----
pdf("./Fig4F-H.pdf", width=12, height=4, onefile = TRUE)
par(mfrow=c(1,3))

for(env_type in c('Tmean', 'DTR', 'STR')){
  if(env_type == 'Tmean'){
    traits_high = all_data %>% filter(Tmean_g == 'High') %>% select(CTmin, CTmax)
    traits_mid = all_data %>% filter(Tmean_g == 'Mid') %>% select(CTmin, CTmax)
    traits_low = all_data %>% filter(Tmean_g == 'Low') %>% select(CTmin, CTmax)
  }else if(env_type == 'DTR'){
    traits_high = all_data %>% filter(DTR_g == 'High') %>% select(CTmin, CTmax)
    traits_mid = all_data %>% filter(DTR_g == 'Mid') %>% select(CTmin, CTmax)
    traits_low = all_data %>% filter(DTR_g == 'Low') %>% select(CTmin, CTmax)
  }else if(env_type == 'STR'){
    traits_high = all_data %>% filter(STR_g == 'High') %>% select(CTmin, CTmax)
    traits_mid = all_data %>% filter(STR_g == 'Mid') %>% select(CTmin, CTmax)
    traits_low = all_data %>% filter(STR_g == 'Low') %>% select(CTmin, CTmax)
  }
  # High gtoup
  H_high <- Hpi(x=traits_high)      # optimal bandwidth estimation
  est_high<- kde(x=traits_high, H=H_high, compute.cont=TRUE)     # kernel density estimation
  # set contour probabilities for drawing contour levels
  cl_high<-contourLevels(est_high, prob=c(0.5, 0.05, 0.001), approx=TRUE)
  
  # Mid gtoup
  H_mid <- Hpi(x=traits_mid)      # optimal bandwidth estimation
  est_mid<- kde(x=traits_mid, H=H_mid, compute.cont=TRUE)     # kernel density estimation
  # set contour probabilities for drawing contour levels
  cl_mid<-contourLevels(est_mid, prob=c(0.5, 0.05, 0.001), approx=TRUE)
  
  # Low gtoup
  H_low <- Hpi(x=traits_low)      # optimal bandwidth estimation
  est_low<- kde(x=traits_low, H=H_low, compute.cont=TRUE)     # kernel density estimation
  # set contour probabilities for drawing contour levels
  cl_low<-contourLevels(est_low, prob=c(0.5, 0.05, 0.001), approx=TRUE)
  
  plot(traits_low, xlim=c(-10, 20), ylim=c(25, 55), 
       ylab="", xlab="", type='n', axes=F,
       main = paste(env_type))
  box(col = 'black') 
  axis(side=1, at=seq(-10,20,5),cex.axis=1.3) # x-axis
  mtext(side=1, line=3, "Critical thermal minimum (\u00B0C)", col="black", font=1,cex=1.1)
  axis(side=2, at=seq(25,55,5),cex.axis=1.3,las=2) # y-axis
  mtext(side=2, line=3, "Critical thermal maximum (\u00B0C)", col="black", font=1,cex=1.1)
  plot(est_high,abs.cont=cl_high[2], labels=c(0.95),labcex=0.75, add=TRUE, lwd=1.2, col=scales::alpha("#CB1B45",0.7))
  plot(est_mid,abs.cont=cl_mid[2], labels=c(0.95),labcex=0.75, add=TRUE, lwd=1.2, col=scales::alpha("#937DC2",0.7))
  plot(est_low,abs.cont=cl_low[2], labels=c(0.95),labcex=0.75, add=TRUE, lwd=1.2, col=scales::alpha('#005CAF',0.7))
  
  legend('bottomright',c(' High',' Mid',' Low'),lty=c(1,1,1),col=scales::alpha(c("#CB1B45",'#937DC2','#005CAF'),0.7),bg='white',bty='n',cex=1.3)
}

dev.off()
