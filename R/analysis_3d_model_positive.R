rm(list=ls())
library(tidyverse)
library(sf)
library(sdmTMB)
library(spdep)
rm(list=ls())
scriptPath <- rstudioapi::getSourceEditorContext()$path 
scriptDir <- dirname(scriptPath)
setwd(file.path(scriptDir, '..'))
source("R/HighstatLibV11.R")

# load data
xdat=read.csv('data/model_3d_input.csv')

# some explorations
ggplot(data=xdat, aes(x=xlon, y=xlat))+
  geom_point(aes(color=log(Abun.rec)))+
  facet_wrap(~year)+
  scale_color_viridis_c()+
  ggtitle('Recruits abundance')

ggplot(data=xdat, aes(x=xlon, y=xlat))+
  geom_point(aes(color=log(Abun.adu)))+
  facet_wrap(~year)+
  scale_color_viridis_c()+
  ggtitle('Adults abundance')

ggplot(data=xdat, aes(x=xlon, y=xlat))+
  geom_point(aes(color=(pp_autumn)))+
  facet_wrap(~year)+
  scale_color_viridis_c()+
  ggtitle('Adults abundance')

ggplot(data=xdat, aes(x=xlon, y=xlat))+
  geom_point(aes(color=log(bottomT_winter)))+
  facet_wrap(~year)+
  scale_color_viridis_c()+
  ggtitle('Adults abundance')



# spatial transformation
xdatsf=xdat%>%
  st_as_sf(coords=c('xlon','xlat'))%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)
xcoords=st_coordinates(xdatsf)
xdat$Xkm=xcoords[,1]/1000
xdat$Ykm=xcoords[,2]/1000

# check distributions
transf.fun=function(xdat, transf.type){
  if(transf.type=='z.score'){
    (xdat - mean(xdat))/sd(xdat)
  }else if (transf.type=='min.max'){
    (xdat - min(xdat))/(max(xdat)-min(xdat))
  }
}

hist(transf.fun(xdat$pp_spring, 'z.score'))
hist(transf.fun(xdat$bottomT_autumn, 'z.score'))
hist(transf.fun(xdat$meandepth, 'z.score'))



# look for spatial autocorrelation
xdatsf2=xdat%>%
  dplyr::group_by(xlon,xlat,year)%>%
  dplyr::summarise(rec=mean(Abun.rec))%>%
  st_as_sf(coords=c('xlon','xlat'))%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)
xdatsf2$logrec=log(xdatsf2$rec+1)

#land.2=st_crop(land,st_buffer(xdatsf2,50000))
#p.full=ggplot()+
#  geom_sf(data=land.2, fill='darkgrey')+
#  geom_sf(data=st_buffer(xdatsf2, 15000), aes(fill=logrec))+
#  scale_fill_gradient2()+
#  labs(fill='log recruitment density')+
#  theme(legend.position = 'bottom')+
#  facet_wrap(~year);p.full
#ggsave('plots/recruit_distribution.png', width = 30, height = 50, units='cm')

# explore spatial autocorr
par(mfrow=c(1,1))
nb <- dnearneigh(st_centroid(xdatsf2), d1=0, d2=5000)
nb_lw <- nb2listw(nb, zero.policy=T)
moran.plot(xdatsf2$logrec, listw=nb_lw)
moran.mc(xdatsf2$logrec, nb_lw, nsim=999)
moran.test(xdatsf2$logrec, nb_lw)
m.t.res=moran.test(xdatsf2$logrec, nb_lw)
table.moran=data.frame(t(m.t.res[["estimate"]]))
table.moran$p.val=m.t.res[["p.value"]]
#write.csv(table.moran, 'Results/release/mira_table.csv', row.names = F)

# when you test the correlation on transformed values, it increases. You may not be really confident in NOT rejecting the null hypothesis of random dispersion. This makes sense, considering how much spatial characteristics are importan for your variable. Find a strategy to not do everything from scratch but allowing to include mesh underneath.

## Models ####
# create mesh 
Loc <- as.matrix(xdat[xdat$Abun.rec>0,c("Xkm", "Ykm")])
D <- dist(Loc)
mesh <- make_mesh(xdat[xdat$Abun.rec>0,], xy_cols = c("Xkm", "Ykm"), cutoff = 10)
plot(mesh)

## dataset for model
x.dat.m=data.frame(pa=ifelse(xdat$Abun.rec>0,1,0))
x.dat.m$Xkm=xdat$Xkm
x.dat.m$Ykm=xdat$Ykm
x.dat.m$rec=xdat$Abun.rec
x.dat.m$logrec=sqrt(xdat$Abun.rec)
x.dat.m$fyear <- as.factor(xdat$year)
x.dat.m$year <- (xdat$year)
x.dat.m$fmonth <- as.factor(xdat$month)
x.dat.m$dep=transf.fun(xdat$meandepth, 'z.score')
x.dat.m$gra=transf.fun(xdat$grain, 'z.score')
x.dat.m$tem.fal=transf.fun(xdat$bottomT_autumn, 'z.score')
x.dat.m$tem.sum=transf.fun(xdat$bottomT_summer, 'z.score')
x.dat.m$tem.spr=transf.fun(xdat$bottomT_spring, 'z.score')
x.dat.m$sal.fal=transf.fun(xdat$sal_autumn, 'z.score')
x.dat.m$sal.sum=transf.fun(xdat$sal_summer, 'z.score')
x.dat.m$sal.spr=transf.fun(xdat$sal_spring, 'z.score')
x.dat.m$ppr.fal=transf.fun(xdat$pp_autumn, 'z.score')
x.dat.m$ppr.sum=transf.fun(xdat$pp_summer, 'z.score')
x.dat.m$ppr.spr=transf.fun(xdat$pp_spring, 'z.score')
x.dat.m=x.dat.m[x.dat.m$pa==1,]

plot(density(x.dat.m$logrec))
plot(density(log(x.dat.m$rec)))
write.csv(x.dat.m, 'data/input_model_pos.csv', row.names = F)

## Presence/Absence
model.metrics=function(x.formula0, xdat, x.spatial, x.spatiotemporal, n.fold, 
                       tv=F, x.formula.tv=NULL){
  
  if(tv==F){
      x.model <- sdmTMB(
  formula = x.formula0,
  data = xdat,
  mesh = mesh,
  family = lognormal(link = "log"),
  spatial=x.spatial,
  time = "year",
  spatiotemporal=x.spatiotemporal
  )
  }else{
    
    x.model <- sdmTMB(
      formula = x.formula0,
      data = xdat,
      mesh = mesh,
      time_varying =x.formula.tv,
      family = lognormal(link = "log"),
      spatial=x.spatial,
      time = "year",
      spatiotemporal=x.spatiotemporal
    ) 
  }
  

# overall check
  x.san=as.data.frame(sanity(x.model))
  x.san=length(x.san[x.san==F])
  
  if(x.san>0){
    
    x.res=data.frame(fail.sanity=x.san, aic=NA,
                     sign.year=NA,
                     sign.fix=NA,
                     cross.val.ll=NA, 
                     rmse=NA,
                     r2=NA,
                     mase=NA)
    return(x.res)
    break
  }
  
  # fixed effects importance
  x.tidy=as.data.frame(tidy(x.model, conf.int=T))
  x.tidy$sign=sign(x.tidy$conf.low)==sign(x.tidy$conf.high)
  x.tidy.fix=x.tidy[-grep('year', x.tidy$term),]
  x.tidy.y=x.tidy[grep('year', x.tidy$term),]
  
  # basic metrics
  x.aic=AIC(x.model)
  x.ll=as.numeric(logLik(x.model))
  
  # predictive skills cross val
  set.seed(123)
  kdat=xdat[sample(1:nrow(xdat)), ] # Shuffle the dataset randomly.
  kdat$folds <- cut(seq(1,nrow(xdat)),breaks=n.fold,labels=FALSE)
  store.kfold=NULL
  for(k in 1:n.fold){
    cat(k)
    # train test
    testIndexes <- which(kdat$folds==k,arr.ind=TRUE)
    testData <- kdat[testIndexes, ]
    trainData <- kdat[-testIndexes, ]
    # new mesh
    kmesh=make_mesh(trainData, xy_cols = c("Xkm", "Ykm"), cutoff = 10)
    # fit model
    x.model.k <- sdmTMB(
      formula = x.formula0,
      data = trainData,
      mesh = kmesh,
      family = lognormal(link = "log"),
      spatial=x.spatial,
      time = "year",
      spatiotemporal=x.spatiotemporal
    )
    
    # predict new data
    p_bin_sim <- predict(x.model.k, newdata = testData)
    testData$pred=x.model.k$family$linkinv(p_bin_sim$est) # exp(p_bin_sim$est)/(1+exp(p_bin_sim$est))
   # testData$pred.pa=ifelse(testData$pred>0.5,1,0)
    
    rmse=sqrt(mean((testData$logrec -testData$pred)^2))
    r2=cor(testData$logrec,testData$pred)^2
    ej=mean(abs(testData$logrec -testData$pred))
    MAD=mean(abs(testData$logrec-mean(testData$logrec)))
    mase=ej/MAD ### MASE
    
    ll=as.numeric(logLik(x.model.k))
    pm=c(rmse, r2, mase, ll)
    store.kfold=rbind(store.kfold, pm)
  }
  cm.k=apply(store.kfold[,1:3], 2, mean)
  ll.k=sum(store.kfold[,4])
  
  # assemble result
  x.res=data.frame(fail.sanity=x.san, aic=x.aic,
                   sign.year=nrow(x.tidy.y[x.tidy.y$sign==T,]),
                   sign.fix=nrow(x.tidy.fix[x.tidy.fix$sign==T,]),
                   cross.val.ll=ll.k, 
                   rmse=cm.k[1],
                   r2=cm.k[2],
                   mase=cm.k[3])
  return(x.res)
}

# positive model the probability to observe recruits. This depends on (i) habitat suitability - depth, salinity, temperature, grainsize - and (ii) timing - how recruits move out from nurseries - time of the year
names(x.dat.m)
# one covariate a time, no random field
x.f.1=logrec ~ 0+fyear +tem.sum
x.f.2=logrec ~ 0+fyear +tem.spr
x.f.3=logrec ~ 0+fyear +tem.fal #
x.f.4=logrec ~ 0+fyear +sal.sum #
x.f.5=logrec ~ 0+fyear +sal.spr
x.f.6=logrec ~ 0+fyear +sal.fal

x.f.7=logrec ~ 0+fyear +ppr.sum #
x.f.8=logrec ~ 0+fyear +ppr.spr
x.f.9=logrec ~ 0+fyear +ppr.fal
x.f.10=logrec ~ 0+fyear +dep #
x.f.11=logrec ~ 0+fyear +s(dep)
x.f.12=logrec ~ 0+fyear +s(sal.fal)
x.elements=ls()
x.elements=x.elements[grep('x.f.',x.elements)]
x.model.list=lapply(x.elements, get)
#names(x.list)=x.elements

model.performance=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                x.spatial = 'off',
                x.spatiotemporal = 'off', n.fold=5)
  x.metrics$id=x.f
  model.performance=rbind(model.performance, x.metrics)  
} 
mp0=model.performance ## retain: temp.spr; sal.spr; sal.fal; ppr.sum; dep
mp0%>%arrange(id)
# test combination of covariates
x.ff.1=logrec ~ 0+fyear +tem.fal + sal.sum #ok
x.ff.2=logrec ~ 0+fyear +tem.fal + sal.fal # ok
x.ff.3=logrec ~ 0+fyear +tem.fal + ppr.sum #ok
x.ff.4=logrec ~ 0+fyear +tem.fal + dep 
x.ff.5=logrec ~ 0+fyear +sal.sum + sal.fal
x.ff.6=logrec ~ 0+fyear +sal.sum + ppr.sum 
x.ff.7=logrec ~ 0+fyear +sal.sum + dep #ok
x.ff.8=logrec ~ 0+fyear +sal.fal + ppr.sum 
x.ff.9=logrec ~ 0+fyear +sal.fal + dep #ok
x.ff.10=logrec ~ 0+fyear +ppr.sum + dep 

x.elements=ls()
x.elements=x.elements[grep('x.ff',x.elements)]
x.model.list=lapply(x.elements, get)

model.performance=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'off',
                          x.spatiotemporal = 'off', n.fold=5)
  x.metrics$id=x.f
  model.performance=rbind(model.performance, x.metrics)  
} 
mp1=model.performance ## the months can be discarded, depth and sal are very important
#mp1$id=seq(1:nrow(mp1))
mp1%>%arrange(rmse)

# keep dep ; ppr.sum ; sal.fal; tem spr

x.fff.1=logrec ~ 0+fyear +s(dep) + sal.sum + ppr.sum
x.fff.2=logrec ~ 0+fyear +s(dep) + sal.sum + tem.fal #
x.fff.3=logrec ~ 0+fyear +sal.sum + tem.fal + ppr.sum
x.fff.4=logrec ~ 0+fyear +sal.sum + tem.fal + ppr.sum + s(dep)

x.elements=ls()
x.elements=x.elements[grep('x.fff',x.elements)]
x.model.list=lapply(x.elements, get)

model.performance=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'off',
                          x.spatiotemporal = 'off', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='off'
  model.performance=rbind(model.performance, x.metrics)  
} 
mp2=model.performance ## the months can be discarded, depth and sal are very important
mp2$id=seq(1:nrow(mp2))
mp2%>%arrange(rmse)

## spatial only on best models
#x.ss.1=logrec ~ 0 +sal.sum + tem.fal + ppr.sum + dep
#x.ss.2=logrec ~ 0 +fyear +dep + sal.sum + tem.fal
#rm(x.ss.2)


#x.elements=ls()
#x.elements=x.elements[grep('x.ss',x.elements)]
#x.model.list=lapply(x.elements, get)

model.performance.s=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'on',
                          x.spatiotemporal = 'off', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='spatial'
  model.performance=rbind(model.performance, x.metrics)  
} 
model.performance.s ## it seems that the simplest model performs bette

#model.performance.s=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'on',
                          x.spatiotemporal = 'iid', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='st_iid'
  model.performance=rbind(model.performance, x.metrics)  
} 

for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'on',
                          x.spatiotemporal = 'rw', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='st_rw'
  model.performance=rbind(model.performance, x.metrics)  
} 

for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'on',
                          x.spatiotemporal = 'ar1', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='st_ar1'
  model.performance=rbind(model.performance, x.metrics)  
} 


model.performance ## it seems that the simplest model performs bette

hist(x.dat.m$logrec)

write.csv(model.performance, 'results/analysis_3d/pos_random_performance.csv', row.names = F)

## more diagnostics on the best model
best.model.pos=sdmTMB(
  formula = logrec ~ 0+fyear +(sal.sum) + tem.fal + s(dep),
  data = x.dat.m,
  mesh = mesh,
  family = lognormal(link = "log"),
  spatial='on',
  time = "year",
  spatiotemporal='off'
)
saveRDS(best.model.pos, 'results/analysis_3d/pos_model.RDS')
library(visreg)
visreg(best.model.pos, xvar = "dep", scale = "response", plot = T)

df.res=data.frame(tidy(best.model.pos, conf.int = TRUE))
p=ggplot(data=df.res[df.res$term!='(Intercept)',])+
  geom_vline(xintercept = 0)+
  geom_point(aes(y=term, x=estimate))+
  geom_errorbar(aes(y=term, xmin=conf.low, xmax=conf.high))
xpreds=simulate(best.model.pos, nsim = 999)
x.dat.m$pred=apply(xpreds, 1, median)


x.dat.m$resid=(apply(xpreds, 1, median))-sqrt(x.dat.m$rec)
ggplot(data=x.dat.m,aes(x=Xkm, y=Ykm))+
  geom_point(aes( color=resid), size=3)+
  facet_wrap(~year)

plot(x.dat.m$pred, x.dat.m$logrec)
abline(a=0,b=1)
par(mfrow=c(2,1))
hist(x.dat.m$pred)
#hist(x.dat.m$logrec)
#hist(sqrt(sqrt(x.dat.m$rec)))
hist(sqrt((x.dat.m$rec)))

ggplot(data=x.dat.m)+
  geom_point(aes(y=resid, x=tem.fal))+
  geom_hline(yintercept=0)


basegrid <- readRDS("results/analysis_3d/grid_pa.RDS")
basegrid$ppr.sum=(basegrid$pp_summer-mean(xdat$pp_summer))/sd(xdat$pp_summer)
coords.grid=as.data.frame(st_coordinates(st_centroid(basegrid)))

pred.grid=as.data.frame(basegrid)%>%
  dplyr::select(year, fyear,dep, sal.fal, ppr.sum, FID)
pred.grid$Xkm=coords.grid$X/1000
pred.grid$Ykm=coords.grid$Y/1000

x.pred.grid=predict(best.model.pos, newdata=pred.grid)
pred.grid$pred.pos=best.model.pos$family$linkinv(x.pred.grid$est)
pred.grid=pred.grid[,c('year', 'FID','pred.pos')]

basegrid=left_join(basegrid, pred.grid, by=c('year','FID'))
basegrid$index=(basegrid$pred.pos)*basegrid$pred
basegrid$index.tv=(basegrid$pred.pos)*basegrid$pred.tv

saveRDS(basegrid, 'results/analysis_3d/grid_pos.RDS')


ggplot(data=basegrid)+
  geom_sf(aes(fill=(pred.pos)*pred), color=NA)+
  facet_wrap(~year)+
  labs(fill='Predicted log Recruit/km2', size='Observed log Recruit/km2')+
  scale_fill_viridis_c()+
  #geom_point(data=x.dat.m, aes(x=Xkm*1000, y=Ykm*1000, size=logrec),color='orange', fill=NA, alpha=0.2)+
  theme(legend.position = 'bottom')

ggsave('results/analysis_3d/plot_prediction.png', width = 50, height = 50, units='cm')

x.dat.m
names(x.dat.m)

# test time varying: it miserably fails
dep + sal.fal + dep
0+fyear +(sal.sum) + tem.fal + ppr.sum + dep

xf1=logrec ~ 0 + sal.sum+ s(dep)
xf1tv= ~ 0 +  tem.fal

mm.tv1=model.metrics(xdat=x.dat.m, x.formula0 = xf1, 
              x.spatial = 'on',
              x.spatiotemporal = 'iid', n.fold=5, tv=T,
              x.formula.tv = xf1tv)

xf2=logrec ~ 0 +  s(dep) + tem.fal
xf2tv= ~ 0 +  sal.sum

mm.tv2=model.metrics(xdat=x.dat.m, x.formula0 = xf2, 
                     x.spatial = 'on',
                     x.spatiotemporal = 'iid', n.fold=5, tv=T,
                     x.formula.tv = xf2tv)

xf3=logrec ~ 0 + sal.sum +tem.fal
xf3tv= ~ 0 +  dep

mm.tv3=model.metrics(xdat=x.dat.m, x.formula0 = xf3, 
                     x.spatial = 'on',
                     x.spatiotemporal = 'iid', n.fold=5, tv=T,
                     x.formula.tv = xf3tv)

xf4=logrec ~ 0 + sal.sum+ s(dep) + ppr.sum
xf4tv= ~ 0 +  tem.fal

mm.tv4=model.metrics(xdat=x.dat.m, x.formula0 = xf4, 
                     x.spatial = 'on',
                     x.spatiotemporal = 'iid', n.fold=5, tv=T,
                     x.formula.tv = xf1tv)


write.csv(rbind(mm.tv1,mm.tv2,mm.tv3, mm.tv4 ), 'results/analysis_3d/pos_random_tv_performance.csv', row.names = F)

best.tv=sdmTMB(
  formula = logrec ~ 0 + sal.sum  + s(dep),
  time_varying = ~ 0 +  tem.fal,
  #time_varying_type = 'rw',
  data = x.dat.m,
  mesh = mesh,
  family = lognormal(link = "log"),
  spatial='on',
  time = "year",
  spatiotemporal='iid'
)
saveRDS(best.tv, 'results/analysis_3d/pos_tv_model.RDS')
summary(best.tv)


nd <- expand.grid(
  tem.fal = seq(min(x.dat.m$tem.fal) + 0.2,
            max(x.dat.m$tem.fal) - 0.2,
            length.out = 50
  ),
  sal.sum=mean(x.dat.m$sal.sum),
  ppr.sum=mean(x.dat.m$ppr.sum),
  dep=mean(x.dat.m$dep),
  year = unique(x.dat.m$year)
)
#nd$fyear=as.factor(nd$year)


p <- predict(best.tv, newdata = nd, se_fit = TRUE, re_form = NA)
p$r.s=ifelse(p$year<2013,'regime1','regime2')
ggplot(p, aes(tem.fal, exp(est),
              # ymin = exp(est - 1.96 * est_se),
              # ymax = exp(est + 1.96 * est_se),
              group = as.factor(year))) +
  geom_line(aes(colour = year, linetype=factor(r.s)), lwd = 1) +
  #geom_ribbon(aes(fill = year), alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) 
  
  
  
  



