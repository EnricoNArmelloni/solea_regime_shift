# This script automatically fit ecological niche models on solemon data. The actual version performs feature selections based on maxent model. On the selected features are fitted a new max ent and a random forest models.
setwd("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift")
rm(list=ls())
library(tidyverse)
library(sf)
library(ncdf4)
'%ni%'=Negate('%in%')
source("C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/R_utilities/netCDF_processing_supp.R")

# sampling stations data ####
occurence <- "C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task3_1/data/solemon_by_haul.csv"
occ <- read.table(occurence, header=TRUE, sep=',')
ch=occ[occ$year==2012,]%>%arrange(haul_number)
haul.points=occ%>%
  distinct(year,month,haul_number,lon,lat,meandepth,sqkm)
xdat_sf=haul.points%>%
  st_as_sf(coords=c('lon','lat'))%>%
  st_set_crs(4326)%>%
  st_transform(3003)

xdat_buffer=st_buffer(xdat_sf, 2500 ) # cerchi con 5km di diametro
xdat_buffer$id=seq(1:nrow(xdat_buffer))
xdat_buffer$xlat=haul.points$lat
xdat_buffer$xlon=haul.points$lon  

ggplot()+
  geom_sf(data=xdat_buffer)



# static layers ####
grainsize=readRDS("C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/grainsize/new_grain.RDS")%>%
  st_set_crs(4326)%>%
  st_transform(3003)%>%
  st_buffer(4631/2)
grainsize$y.id=seq(1:nrow(grainsize))
xint=data.frame(st_intersects(grainsize, xdat_buffer))
names(xint)[1]='y.id'
xint=xint%>%left_join(grainsize, by='y.id')%>%
  dplyr::group_by(col.id)%>%
  dplyr::summarise(grain=mean(grain_17))
names(xint)[1]='id'
xdat_buffer=xdat_buffer%>%
  left_join(xint, by=c('id'))


## spatiotemporal covariates ####
gsa=read_sf("C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/gsa/GSAs_simplified")%>%
  st_set_crs(4326)%>%
  dplyr::filter(SMU_CODE==17)

temp.dat=open_nc(yrs=2007:2019,
        xvar='bottomT',
        stratification='no' ,
        myfolder='bottomT',
        region='Adriatic',
        time_res='monthly',
        basedir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/oceanography_physic",
        croplayer = gsa)

temp.dat$month=as.numeric(temp.dat$month)
temp.dat$season=ifelse(temp.dat$month %in% 1:4, 'winter',
                       ifelse(temp.dat$month %in% 5:6, 'spring',
                              ifelse(temp.dat$month %in% 7:9, 'summer','autumn')))
bottomT=temp.dat%>%
  dplyr::filter(!is.na(xvar))%>%
  dplyr::group_by(lon,lat,year,season)%>%
  dplyr::summarise(bottomT=mean(xvar))
bottomT=bottomT%>%
  st_as_sf(coords=c('lon','lat'))%>%
  st_set_crs(4326)%>%
  st_transform(3003)%>%
  st_buffer(4631/2)


sal.dat=open_nc(yrs=2007:2019,
                 xvar='sal',
                 stratification='bottom' ,
                 myfolder='salinity',
                 region='Adriatic',
                 time_res='monthly',
                 basedir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/oceanography_physic",
                 croplayer = gsa)
sal.dat$month=as.numeric(sal.dat$month)
sal.dat$season=ifelse(sal.dat$month %in% 1:4, 'winter',
                       ifelse(sal.dat$month %in% 5:6, 'spring',
                              ifelse(sal.dat$month %in% 7:9, 'summer','autumn')))
bottomSAL=sal.dat%>%
  dplyr::filter(!is.na(xvar))%>%
  dplyr::group_by(lon,lat,year,season)%>%
  dplyr::summarise(sal=mean(xvar))
bottomSAL=bottomSAL%>%
  st_as_sf(coords=c('lon','lat'))%>%
  st_set_crs(4326)%>%
  st_transform(3003)%>%
  st_buffer(4631/2)


pp.dat=open_nc(yrs=2007:2019,
        xvar='prim_prod',
        stratification='bottom' ,
        myfolder='prim_prod',
        region='Adriatic',
        time_res='monthly',
        basedir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/oceanography_biogeo",
        croplayer = gsa)
pp.dat$month=as.numeric(pp.dat$month)
pp.dat$season=ifelse(pp.dat$month %in% 1:4, 'winter',
                      ifelse(pp.dat$month %in% 5:6, 'spring',
                             ifelse(pp.dat$month %in% 7:9, 'summer','autumn')))
bottomPP=pp.dat%>%
  dplyr::filter(!is.na(xvar))%>%
  dplyr::group_by(lon,lat,year,season)%>%
  dplyr::summarise(pp=mean(xvar))
bottomPP=bottomPP%>%
  st_as_sf(coords=c('lon','lat'))%>%
  st_set_crs(4326)%>%
  st_transform(3003)%>%
  st_buffer(4631/2)


xyrs=2007:2019#as.numeric(unique(idat$year))
xyrs=sort(xyrs)
est_store=list()

for(j in 1:length(xyrs)){
    
  cat(xyrs[j])
    # load species data and subset on solemon correct strata
    xdat.haul=xdat_buffer[xdat_buffer$year==xyrs[j], ]
    xdat.haul=xdat.haul[xdat.haul$month>6,]
    xdat.haul$x.id=seq(1:nrow(xdat.haul))
    
    # bottom T
    xdat.bottom=bottomT[bottomT$year==xyrs[j],]
    xdat.bottom$y.id=seq(1:nrow(xdat.bottom))
    xint=data.frame(st_intersects(xdat.bottom, xdat.haul))
    names(xint)[1]='y.id'
    xint=xint%>%left_join(xdat.bottom, by='y.id')%>%
      dplyr::group_by(col.id, season)%>%
      dplyr::summarise(bottomT=mean(bottomT))
    names(xint)[1]='x.id'
    xint=xint%>%pivot_wider(names_from = season, values_from = bottomT)
    names(xint)[2:5]=paste0('bottomT_',names(xint[2:5]))
    xdat.haul=xdat.haul%>%
      left_join(xint, by=c('x.id'))
    
    # sal
    xdat.sal=bottomSAL[bottomSAL$year==xyrs[j],]
    xdat.sal$y.id=seq(1:nrow(xdat.sal))
    xint=data.frame(st_intersects(xdat.sal, xdat.haul))
    names(xint)[1]='y.id'
    xint=xint%>%left_join(xdat.sal, by='y.id')%>%
      dplyr::group_by(col.id, season)%>%
      dplyr::summarise(sal=mean(sal))
    names(xint)[1]='x.id'
    xint=xint%>%pivot_wider(names_from = season, values_from = sal)
    names(xint)[2:5]=paste0('sal_',names(xint[2:5]))
    xdat.haul=xdat.haul%>%
      left_join(xint, by=c('x.id'))
    
    # bottom T
    xdat.pp=bottomPP[bottomPP$year==xyrs[j],]
    xdat.pp$y.id=seq(1:nrow(xdat.pp))
    xint=data.frame(st_intersects(xdat.pp, xdat.haul))
    names(xint)[1]='y.id'
    xint=xint%>%left_join(xdat.pp, by='y.id')%>%
      dplyr::group_by(col.id, season)%>%
      dplyr::summarise(pp=mean(pp))
    names(xint)[1]='x.id'
    xint=xint%>%pivot_wider(names_from = season, values_from = pp)
    names(xint)[2:5]=paste0('pp_',names(xint[2:5]))
    xdat.haul=xdat.haul%>%
      left_join(xint, by=c('x.id'))
    
    xdat.haul=as.data.frame(xdat.haul)%>%
      dplyr::select(-geometry)
    
    est_store[[j]]=xdat.haul
    
}
hauls_covariates=plyr::ldply(est_store)    
hauls_covariates=hauls_covariates%>%dplyr::select(-id,-x.id)    

# paste solea data
setwd("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift")
library(readxl)
rec.data=read_excel("data/station_index_below20cm.xlsx")
rec.data$year=as.numeric(str_remove(str_remove(str_remove(str_remove(rec.data$Survey,'SOLEMON'),'_b'),'OTT'),'NOVEMBRE'))
rec.data$StationDate=as.Date(rec.data$StationDate,origin="1900-01-01")
rec.data$month=lubridate::month(rec.data$StationDate)
rec.data=rec.data[rec.data$month>3,]
rec.data=rec.data[,c('year','Station','AbunIndex')]
names(rec.data)[c(2,3)]=c('haul_number', 'Abun.rec')
rec.data=rec.data[rec.data$year!=2009,]

rec.data.2009=read_excel("data/station_index_below20cm_y2009.xlsx")
rec.data.2009$year=as.numeric(str_remove(str_remove(str_remove(str_remove(rec.data.2009$Survey,'SOLEMON'),'_b'),'OTT'),'NOVEMBRE'))
rec.data.2009$StationDate=as.Date(rec.data.2009$StationDate,origin="1900-01-01")
rec.data.2009$month=lubridate::month(rec.data.2009$StationDate)
rec.data.2009=rec.data.2009[rec.data.2009$month>3,]
rec.data.2009=rec.data.2009[,c('year','Station','AbunIndex')]
names(rec.data.2009)[c(2,3)]=c('haul_number', 'Abun.rec')
rec.data.2009=rec.data.2009[rec.data.2009$year==2009,]

rec.data=rbind(rec.data,rec.data.2009)%>%
  arrange(year, haul_number)



adu.data=read_excel("data/station_index_above20cm.xlsx")
adu.data$year=as.numeric(str_remove(str_remove(str_remove(str_remove(adu.data$Survey,'SOLEMON'),'_b'),'OTT'),'NOVEMBRE'))
adu.data$StationDate=as.Date(adu.data$StationDate,origin="1900-01-01")
adu.data$month=lubridate::month(adu.data$StationDate)
adu.data=adu.data[rec.data$month>3,]
adu.data=adu.data[,c('year','Station','AbunIndex', 'BiomIndex')]
names(adu.data)[c(2:4)]=c('haul_number', 'Abun.adu','Biom.adu')
adu.data=adu.data[adu.data$year!=2009,]


adu.data.2009=read_excel("data/station_index_above20cm_y2009.xlsx")
adu.data.2009$year=as.numeric(str_remove(str_remove(str_remove(str_remove(adu.data.2009$Survey,'SOLEMON'),'_b'),'OTT'),'NOVEMBRE'))
adu.data.2009$StationDate=as.Date(adu.data.2009$StationDate,origin="1900-01-01")
adu.data.2009$month=lubridate::month(adu.data.2009$StationDate)
adu.data.2009=adu.data.2009[adu.data.2009$month>3,]
adu.data.2009=adu.data.2009[,c('year','Station','AbunIndex', 'BiomIndex')]
names(adu.data.2009)[c(2:4)]=c('haul_number', 'Abun.adu','Biom.adu')
adu.data.2009=adu.data.2009[adu.data.2009$year==2009,]

adu.data=rbind(adu.data,adu.data.2009)%>%
  arrange(year, haul_number)


solemon.solea=full_join(rec.data, adu.data,by = join_by(year, haul_number))

solemon.solea=solemon.solea[solemon.solea$year %in% xyrs,]

# merge all
solea.dat=right_join(hauls_covariates, solemon.solea, by=c('year','haul_number'))
solea.dat=solea.dat%>%na.omit()
write.csv(solea.dat, 'data/model_input.csv', row.names = F)





