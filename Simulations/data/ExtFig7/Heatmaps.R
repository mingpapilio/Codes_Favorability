##
library(ggplot2)
library(ggpubr)
theme_set(theme_minimal())
library(data.table)
library(ggpmisc)
##
hVol_summary<- fread("hvol_summary.txt")
a<- ggplot(hVol_summary, aes(Tmean, Tsd_short, fill= CTmax_CTmin_hv))+
  geom_tile()+
  scale_fill_viridis_b(limits = c(5, 17.5), n.breaks=10, option = "plasma", name= "")+
  ggtitle("Hypervolume")+
  rremove("ylab") + rremove("xlab")+
  rremove("ylab") + rremove("xlab")+
  theme(plot.title=element_text(size=11, face="italic"),
        legend.key.height= grid::unit(1.5, "cm"),
        legend.key.width= grid::unit(0.3, "cm"))
##
CTrange_summary<- fread("range_summary.txt")
b<- ggplot(CTrange_summary, aes(Tmean, Tsd_short, fill= CTRangeAvg))+
  geom_tile()+
  scale_fill_viridis_b(limits = c(17, 23), n.breaks=10, option = "plasma", name= "")+
  ggtitle("Thermal tolerance range")+
  rremove("ylab") + rremove("xlab")+
  theme(plot.title=element_text(size=11, face="italic"),
        legend.key.height= grid::unit(1.5, "cm"),
        legend.key.width= grid::unit(0.3, "cm"))
nspc_summary<- fread("nspc_summary.txt")
c<- ggplot(nspc_summary, aes(Tmean, Tsd_short, fill= NumSpc))+
  geom_tile()+
  scale_fill_viridis_b(limits = c(1000, 3000), n.breaks=10, option = "plasma", name= "")+
  ggtitle("Number of species")+
  rremove("ylab") + rremove("xlab")+
  theme(plot.title=element_text(size=11, face="italic"),
        legend.key.height= grid::unit(1.5, "cm"),
        legend.key.width= grid::unit(0.3, "cm"))
Span_summary<- fread("span_summary.txt")
d<- ggplot(Span_summary, aes(Tmean, Tsd_short, fill= SpanAvg))+
  geom_tile()+
  scale_fill_viridis_b(limits = c(5e4, 7e4), n.breaks=10, option = "plasma", name= "")+
  ggtitle("Average lifespan")+
  rremove("ylab") + rremove("xlab")+
  theme(plot.title=element_text(size=11, face="italic"),
        legend.key.height= grid::unit(1.5, "cm"),
        legend.key.width= grid::unit(0.3, "cm"))

fig<- ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d"))
fig<- annotate_figure(fig, left=text_grob("Thermal variability (°C)", rot= 90, size= 14), bottom= text_grob("Mean ambient temperature (°C)", size= 14))
ggsave("HeatMap.pdf",fig, width = 10, height= 9, units="in")
 