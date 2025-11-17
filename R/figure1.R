rm(list=ls())
scriptPath <- rstudioapi::getSourceEditorContext()$path 
scriptDir <- dirname(scriptPath)
setwd(file.path(scriptDir, '..'))
library(readxl)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(ggspatial)
library(prettymapr)

# Fig 1 A
gsas=read_sf("C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/gsa/GSAs_simplified")%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)
gsa=gsas[gsas$SMU_CODE %in%c(17),]
xcountry=ne_countries(country = c("italy"), scale = "medium")%>%
  st_as_sf%>%
  st_set_crs(4326)%>%
  st_transform(3003)

spat.pred=readRDS("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/results/spatial_predictions.RDS")



solemon=read_sf('C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/Data/shapefiles/StrataITA_SVN')%>%
  st_set_crs(4326)%>%
  st_transform(3003)

closure6=st_buffer(gsa, -6*1000*1.852)
closure3=st_buffer(gsa, -3*1000*1.852)
ch=st_difference(gsa, closure6)

inset=ggplot()+
  annotation_map_tile() +
  geom_sf(data=gsas, fill=NA)+
  geom_sf(data=gsa, aes(fill=as.character(SMU_CODE)))+
  #geom_sf(data=FMZ, aes(fill=fmz), alpha=0.5)+
  #geom_sf(data=xcountry, fill='darkgrey')+
  annotation_scale() +
  labs(fill='Area')+
  #geom_sf(data=idat.spat)+
  scale_fill_manual(values=c('tomato1'))+
  theme(legend.position = 'bottom');inset

closure3=st_crop(closure3, t1)
closure6=st_crop(closure6, t1)

map1=ggplot()+
  annotation_map_tile() +
  geom_sf(data=gsa, fill=NA)+
  #geom_sf(data=gsa15, fill='tomato1', color='tomato1')+
  #geom_sf(data=FMZ, aes(fill=fmz), alpha=0.2)+
  #geom_sf(data=xcountry, fill='darkgrey')+
  annotation_scale() +
  geom_sf(data=solemon, alpha=0.4)+
  #geom_sf(data=closure3, fill=NA, color='red')+
  #geom_sf(data=closure6, fill=NA, color='red')+
  labs(fill='Area')+
  scale_fill_manual(values=c('seashell','tomato1'));map1

map3=ggpubr::ggarrange(map1,map2, ncol=1, heights = c(0.85,2))
ggsave(plot=map3,'msfd_malta_r2/results/supporting_maps/study_area.jpeg', width = 15, height = 20, units = 'cm', dpi=150)



### Fig 2 d-g
dat.sol=read_csv("data/solea_tot.csv")

p2d=ggplot(data=dat.sol, aes(x=year, y=runoff_ms))+
  geom_smooth(color='black')+
  geom_point()+
  geom_line()+
  theme_bw()+
  ylab(paste('Runoff Po River (',expression(m^3) ,'/second)'))+
  xlab('Year')

p2e=ggplot(data=dat.sol, aes(x=year, y=PP ))+
  geom_smooth(color='black')+
  geom_point()+
  geom_line()+
  theme_bw()+
  ylab(paste('Primary Production (C/unit volume)'))+
  xlab('Year')

p2g=ggplot(data=dat.sol, aes(x=year, y=effort ))+
  geom_smooth(color='black')+
  geom_point()+
  geom_line()+
  theme_bw()+
  ylab(expression('Fishing capacity F'[INT]))+
  xlab('Year')

p2f=ggplot(data=dat.sol, aes(x=year, y=bottom_t ))+
  geom_smooth(color='black')+
  geom_point()+
  geom_line()+
  theme_bw()+
  ylab(paste('Sea Bottom Temperature (C°)'))+
  xlab('Year')

p2df=ggpubr::ggarrange(p2d,p2e,p2f,p2g,ncol=1)



## 2 b c
t1=spat.pred[spat.pred$Regime1>0,]%>%
  dplyr::group_by(x)%>%
  dplyr::summarise(Regime1=mean(Regime1))
t2=spat.pred[spat.pred$Regime2>0,]%>%
  dplyr::group_by(x)%>%
  dplyr::summarise(Regime1=mean(Regime2))

p2b=ggplot()+
  annotation_map_tile() +
  #geom_sf(data=gsa, fill=NA)+
  #geom_sf(data=gsa15, fill='tomato1', color='tomato1')+
  #geom_sf(data=FMZ, aes(fill=fmz), alpha=0.2)+
  #geom_sf(data=xcountry, fill='darkgrey')+
  annotation_scale() +
  geom_sf(data=solemon, alpha=0.4)+
  geom_sf(data=t1,aes(fill=log(Regime1+1)), color=NA)+
  #geom_sf(data=closure3, fill=NA, color='red')+
  geom_sf(data=closure6, fill=NA, color='red')+
  coord_sf(xlim = c(1400000, 1720000 ), ylim = c(5150000 , 5500000 ), expand = FALSE)+
  labs(fill='Area')+scale_fill_viridis_c();p2b

p2c=ggplot()+
  annotation_map_tile() +
  #geom_sf(data=gsa, fill=NA)+
  #geom_sf(data=gsa15, fill='tomato1', color='tomato1')+
  #geom_sf(data=FMZ, aes(fill=fmz), alpha=0.2)+
  #geom_sf(data=xcountry, fill='darkgrey')+
  annotation_scale()+
  geom_sf(data=solemon, alpha=0.4)+
  geom_sf(data=t2,aes(fill=log(Regime1+1)), color=NA)+
  #geom_sf(data=closure3, fill=NA, color='red')+
  geom_sf(data=closure6, fill=NA, color='red')+
  coord_sf(xlim = c(1400000, 1720000 ), ylim = c(5150000 , 5500000 ), expand = FALSE)+
  labs(fill='Area')+scale_fill_viridis_c()

ggpubr::ggarrange(p2b,p2c)


p2b=ggplot()+
  annotation_map_tile() +
  #geom_sf(data=gsa, fill=NA)+
  #geom_sf(data=gsa15, fill='tomato1', color='tomato1')+
  #geom_sf(data=FMZ, aes(fill=fmz), alpha=0.2)+
  #geom_sf(data=xcountry, fill='darkgrey')+
  annotation_scale()+
  geom_sf(data=solemon, alpha=0.4)+
  geom_sf(data=spat.pred,aes(fill=diff), color=NA)+
  #geom_sf(data=closure3, fill=NA, color='red')+
  geom_sf(data=ch, fill='seashell', alpha=0.2)+
  coord_sf(xlim = c(1350000, 1600000 ), ylim = c(5150000 , 5550000 ), expand = FALSE)+
  labs(fill='Increase in juveniles density')+scale_fill_viridis_c()+
  theme(legend.position = 'bottom');p2b

ggsave(plot=map1, 'plots/Figure1/Figure1_a.png', width=20, height=14, units='cm', dpi=300)
ggsave(plot=p2b, 'plots/Figure1/Figure1.png', width=10, height=14, units='cm', dpi=300)
ggsave(plot=p2df, 'plots/Figure1/Figure1_d.png', width=12, height=28, units='cm', dpi=300)
ggsave(plot=inset, 'plots/Figure1/Figure1_ins.png', width=15, height=10, units='cm', dpi=300)

