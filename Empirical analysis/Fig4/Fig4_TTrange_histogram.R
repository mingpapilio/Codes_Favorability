library(dplyr)
library(ggplot2)

df.all <- read.csv('./FAV_hv_data.csv',header = TRUE)
all_data = df.all %>% dplyr::select(TTrange, Location, STmean, DTR, STR, Group) # select traits needed for analysis

#### Fig.4C,4F,4I (histogram) ####
## grouping species into high, mid and low STmean, DTR and STR
## 'high, mid and low' is defined by the average temperature of each assemblage from high to low, which is illustrated in Fig.4B,4E,4H
all_data = all_data %>% mutate(Tmean_g = ifelse(Group == 'M1'|Group == 'M2'|Group == 'T1','High', 
                                                ifelse(Group == 'T2'|Group == 'C1'|Group == 'T3', 'Mid', 'Low')))
all_data = all_data %>% mutate(DTR_g = ifelse(Group == 'C1'|Group == 'C3'|Group == 'T1', 'High',
                                              ifelse(Group == 'C2'|Group == 'C4'|Group == 'T2', 'Mid','Low')))
all_data = all_data %>% mutate(STR_g = ifelse(Group == 'C4'|Group == 'C3'|Group == 'C2', 'High',
                                              ifelse(Group == 'C1'|Group == 'T2'|Group == 'T3', 'Mid','Low')))

#---- Plotting histogram ----
## global option ##
xmax<- 50 # Maximal x value in the plot
ymax<- 120 # Maximal y value in the plot
binwd<- 3.5 # Width of each bin
axis_size<- 12
text_size<- 13
title_size<- 14
theme_set(theme_bw())

###---- Fig.4C (Tmean)
colors <- c("High Tmean" = "#CB1B45", "Low Tmean" = "#005CAF")
data_low_t = all_data %>% filter(Tmean_g == 'Low') %>% select(TTrange)
data_high_t = all_data %>% filter(Tmean_g == 'High') %>% select(TTrange)
assign('Tmean', ggplot()+
         coord_cartesian(xlim=c(20, 50), ylim = c(0, ymax))+
         ggtitle('')+
         xlab("Thermal tolerance range (°C)")+ ylab("Frequency")+
         geom_histogram(aes(x= data_high_t$TTrange, fill="High Tmean"),binwidth= binwd, na.rm = T, alpha=0.5, color='grey80', size=0.1)+
         geom_histogram(aes(x= data_low_t$TTrange, fill="Low Tmean"),binwidth= binwd, na.rm = T, alpha=0.5, color='grey80', size=0.1)+
         theme(axis.title= element_text(size= title_size), 
               axis.text= element_text(size= axis_size),
               legend.text=element_text(size=axis_size),
               legend.position = "inside", 
               legend.position.inside=c(.15,.9))+
         labs(fill = "")+
         scale_fill_manual(values = colors)+ 
         bgcolor("white"))
Tmean

###---- Fig.4F (DTR)
colors <- c("High DTR" = "#CB1B45", "Low DTR" = "#005CAF")
data_low_d = all_data %>% filter(DTR_g == 'Low') %>% select(TTrange)
data_high_d = all_data %>% filter(DTR_g == 'High') %>% select(TTrange)
assign('DTR', ggplot()+
         coord_cartesian(xlim=c(20, 50), ylim = c(0, ymax))+
         ggtitle('')+
         xlab("Thermal tolerance range (°C)")+ ylab("Frequency")+
         geom_histogram(aes(x= data_high_d$TTrange, fill="High DTR"),binwidth= binwd, na.rm = T, alpha=0.5, color='grey80', size=0.1)+
         geom_histogram(aes(x= data_low_d$TTrange, fill="Low DTR"),binwidth= binwd, na.rm = T, alpha=0.5, color='grey80', size=0.1)+
         theme(axis.title= element_text(size= title_size), 
               axis.text= element_text(size= axis_size),
               legend.text=element_text(size=axis_size),
               legend.position = "inside", 
               legend.position.inside=c(.15,.9))+
         labs(fill = "")+
         scale_fill_manual(values = colors)+ 
         bgcolor("white"))
DTR

###---- Fig.4I (STR)
colors <- c("High STR" = "#CB1B45", "Low STR" = "#005CAF")
data_low_s = all_data %>% filter(STR_g == 'Low') %>% select(TTrange)
data_high_s = all_data %>% filter(STR_g == 'High') %>% select(TTrange)
assign('STR', ggplot()+
         coord_cartesian(xlim=c(20, 50), ylim = c(0, ymax))+
         ggtitle('')+
         xlab("Thermal tolerance range (°C)")+ ylab("Frequency")+
         geom_histogram(aes(x= data_high_s$TTrange, fill="High STR"),binwidth= binwd, na.rm = T, alpha=0.5, color='grey80', size=0.1)+
         geom_histogram(aes(x= data_low_s$TTrange, fill="Low STR"),binwidth= binwd, na.rm = T, alpha=0.5, color='grey80', size=0.1)+
         theme(axis.title= element_text(size= title_size), 
               axis.text= element_text(size= axis_size),
               legend.text=element_text(size=axis_size),
               legend.position = "inside", 
               legend.position.inside=c(.15,.9))+
         labs(fill = "")+
         scale_fill_manual(values = colors)+ 
         bgcolor("white"))
STR

blank = ggplot() + theme_void() # Add blank plot so the plots can be composed easier

### combine
plot_ls <- mget(c('Tmean','DTR','STR'))
fig_comp1<- ggarrange(plot_ls[[1]],plot_ls[[2]],plot_ls[[3]],blank,ncol=2,nrow=2)
fig_comp1
ggsave("./TTrange_hist.pdf",fig_comp1, width = 10, height= 10, units="in")
