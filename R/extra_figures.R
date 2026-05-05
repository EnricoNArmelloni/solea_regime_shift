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
#gsas=read_sf("C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/gsa/GSAs_simplified")%>%
gsas=read_sf('../other_data/GSAs_simplified')%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)
gsa=gsas[gsas$SMU_CODE %in%c(17),]

xcountry=ne_countries(continent=c('Europe', 'Africa', 'Asia'), scale = "medium")%>%
  st_as_sf%>%
  st_set_crs(4326)%>%
  st_transform(3003)%>%
  st_crop(st_buffer(gsas,200*1000))%>%
  st_union

sarea=ggplot()+
  geom_sf(data=gsas, fill='grey90')+
  geom_sf(data=gsa, fill='yellow')+
  geom_sf(data=xcountry, fill='grey20')+
  geom_sf_text(data=gsas, aes(label=str_remove(F_GSA_LIB, 'GSA ')))+
  annotation_scale()+
  xlab('')+
  ylab('')+
  theme(panel.background = element_rect(fill='white'))

ggsave(plot=sarea, 'results/figS1.jpeg', width = 25, height = 15, units='cm', dpi=500)

### Fig 2 d-g
dat.sol=read_csv("data/model_2d_input.csv")
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
  ylab(paste('Sea Bottom Temperature (C)'))+
  xlab('Year')

p2df=ggpubr::ggarrange(p2d,p2e,p2f,p2g, labels=c('a)','b)','c)','d)'))
ggsave(plot=p2df,'results/covars.jpeg', width = 18, height = 18, units = 'cm', dpi=500)

