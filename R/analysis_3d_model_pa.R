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
Loc <- as.matrix(xdat[,c("Xkm", "Ykm")])
D <- dist(Loc)
mesh <- make_mesh(xdat, xy_cols = c("Xkm", "Ykm"), cutoff = 10)
plot(mesh)

## dataset for model
x.dat.m=data.frame(pa=ifelse(xdat$Abun.rec>0,1,0))
x.dat.m$Xkm=xdat$Xkm
x.dat.m$Ykm=xdat$Ykm
x.dat.m$logrec=log(xdat$Abun.rec+1)
x.dat.m$fyear <- as.factor(xdat$year)
x.dat.m$year <- (xdat$year)
x.dat.m$fmonth <- as.factor(xdat$month)
x.dat.m$dep=transf.fun(xdat$meandepth, 'z.score')
x.dat.m$gra=transf.fun(xdat$grain, 'z.score')
x.dat.m$tem.fal=transf.fun(xdat$bottomT_autumn, 'z.score')
x.dat.m$sal.fal=transf.fun(xdat$sal_autumn, 'z.score')


## collinearity expl
corvif(x.dat.m[,c('dep', 'gra', 'tem.fal', 'sal.fal')]) # all good
write.csv(x.dat.m, 'data/input_model_pa.csv', row.names = F)

## Presence/Absence
model.metrics=function(x.formula0, xdat, x.spatial, x.spatiotemporal, n.fold, 
                       tv=F, x.formula.tv=NULL){
  
  if(tv==F){
      x.model <- sdmTMB(
  formula = x.formula0,
  data = xdat,
  mesh = mesh,
  family = binomial(link = "logit"),
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
      family = binomial(link = "logit"),
      spatial=x.spatial,
      time = "year",
      spatiotemporal=x.spatiotemporal
    ) 
  }
  

# overall check
  x.san=as.data.frame(sanity(x.model))
  x.san=length(x.san[x.san==F])
  
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
      family = binomial(link = "logit"),
      spatial=x.spatial,
      time = "year",
      spatiotemporal=x.spatiotemporal
    )
    
    # predict new data
    p_bin_sim <- predict(x.model.k, newdata = testData)
    testData$pred=x.model.k$family$linkinv(p_bin_sim$est) # exp(p_bin_sim$est)/(1+exp(p_bin_sim$est))
    testData$pred.pa=ifelse(testData$pred>0.5,1,0)
    
    # performance metric
    cm=caret::confusionMatrix(table(testData$pred.pa, testData$pa))
    cm=as.numeric(c(cm$overall[1],cm$byClass[c(1,2,5)]))
    ll=as.numeric(logLik(x.model.k))
    pm=c(cm, ll)
    store.kfold=rbind(store.kfold, pm)
  }
  cm.k=apply(store.kfold[,1:4], 2, mean)
  ll.k=sum(store.kfold[,5])
  
  # assemble result
  x.res=data.frame(fail.sanity=x.san, aic=x.aic,
                   sign.year=nrow(x.tidy.y[x.tidy.y$sign==T,]),
                   sign.fix=nrow(x.tidy.fix[x.tidy.fix$sign==T,]),
                   cross.val.ll=ll.k, 
                   accuracy=cm.k[1],
                   sensitivity=cm.k[2],
                   specificity=cm.k[3],
                   precision=cm.k[4])
  return(x.res)
}

# PA model the probability to observe recruits. This depends on (i) habitat suitability - depth, salinity, temperature, grainsize - and (ii) timing - how recruits move out from nurseries - time of the year

# one covariate a time, no random field
x.f.1=pa ~ 0+fyear +fmonth
x.f.2=pa ~ 0+fyear +dep
x.f.3=pa ~ 0+fyear +gra
x.f.4=pa ~ 0+fyear +tem.fal
x.f.5=pa ~ 0+fyear +sal.fal

x.elements=ls()
x.elements=x.elements[grep('x.f.',x.elements)]
x.model.list=lapply(x.elements, get)
#names(x.list)=x.elements

model.performance=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                x.spatial = 'off',
                x.spatiotemporal = 'off', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='off'
  model.performance=rbind(model.performance, x.metrics)  
} 
mp0=model.performance ## the months can be discarded, depth and sal are very important

# test combination of covariates
x.ff.1=pa ~ 0+fyear +dep + sal.fal
x.ff.2=pa ~ 0+fyear +dep + gra
x.ff.3=pa ~ 0+fyear +dep + tem.fal # worst
x.ff.4=pa ~ 0+fyear +sal.fal + gra
x.ff.5=pa ~ 0+fyear +sal.fal + tem.fal

x.elements=ls()
x.elements=x.elements[grep('x.ff',x.elements)]
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
mp1=model.performance ## the months can be discarded, depth and sal are very important

x.fff.1=pa ~ 0+fyear +dep + sal.fal + gra
x.fff.2=pa ~ 0+fyear +sal.fal + tem.fal + gra
x.fff.3=pa ~ 0+fyear +dep + tem.fal + gra
x.fff.4=pa ~ 0+fyear +dep + sal.fal+ tem.fal + gra

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

model.performance$id=x.elements
met1=as.numeric(quantile(model.performance$accuracy, 0.5))
met2=as.numeric(quantile(model.performance$sensitivity, 0.5))
met3=as.numeric(quantile(model.performance$specificity, 0.5))
met4=as.numeric(quantile(model.performance$precision, 0.5))

best.models=model.performance[model.performance$accuracy>=met1 &
                    model.performance$sensitivity>=met2 &
                    model.performance$specificity>=met3&
                    model.performance$precision>=met4,]



## spatial only on best models
x.s.1=pa ~ 0+fyear +dep + sal.fal
x.s.2=pa ~ 0+fyear +dep + sal.fal + gra
x.s.3=pa ~ 0+fyear +dep + sal.fal+ tem.fal 


x.elements=ls()
x.elements=x.elements[grep('x.s',x.elements)]
x.model.list=lapply(x.elements, get)

model.performance.s=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'off',
                          x.spatiotemporal = 'off', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='spatial'
  model.performance.s=rbind(model.performance.s, x.metrics)  
} 
model.performance.s ## it seems that the simplest model performs bette

#model.performance.s=NULL
for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'on',
                          x.spatiotemporal = 'iid', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='st_iid'
  model.performance.s=rbind(model.performance.s, x.metrics)  
} 

for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'on',
                          x.spatiotemporal = 'rw', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='st_rw'
  model.performance.s=rbind(model.performance.s, x.metrics)  
} 

for(x.f in 1:length(x.model.list)){
  x.metrics=model.metrics(xdat=x.dat.m, x.formula0 = x.model.list[[x.f]], 
                          x.spatial = 'on',
                          x.spatiotemporal = 'ar1', n.fold=5)
  x.metrics$cov.type=x.f
  x.metrics$random='st_ar1'
  model.performance.s=rbind(model.performance.s, x.metrics)  
} 


model.performance.s ## it seems that the simplest model performs bette
write.csv(model.performance.s, 'results/analysis_3d/pa_fix_performance.csv', row.names = F)


best.models=model.performance.s[model.performance.s$fail.sanity==0,]
best.models$accuracy=round(best.models$accuracy, digits=2)
best.models$sensitivity=round(best.models$sensitivity, digits=2)
best.models$specificity=round(best.models$specificity, digits=2)
best.models$precision=round(best.models$precision, digits=2)

best.models


## Best model: fit and more diagnostics ####
best.model.pa=sdmTMB(
  formula = pa ~ 0 +fyear +dep + sal.fal+ tem.fal,
  data = x.dat.m,
  mesh = mesh,
  family = binomial(link = "logit"),
  spatial='on',
  time = "year",
  spatiotemporal='rw'
)
saveRDS(best.model.pa, 'results/analysis_3d/pa_model.RDS')



# test time varying
xf1=pa ~ 0 + sal.fal+ tem.fal
xf1tv= ~ 0 +  dep

mm.tv1=model.metrics(xdat=x.dat.m, x.formula0 = xf1, 
              x.spatial = 'on',
              x.spatiotemporal = 'rw', n.fold=5, tv=T,
              x.formula.tv = xf1tv)

xf2=pa ~ 0 +  tem.fal + dep
xf2tv= ~ 0 +  sal.fal

mm.tv2=model.metrics(xdat=x.dat.m, x.formula0 = xf2, 
                     x.spatial = 'on',
                     x.spatiotemporal = 'rw', n.fold=5, tv=T,
                     x.formula.tv = xf2tv)

xf3=pa ~ 0 + sal.fal  + dep
xf3tv= ~ 0 +  tem.fal

mm.tv3=model.metrics(xdat=x.dat.m, x.formula0 = xf3, 
                     x.spatial = 'on',
                     x.spatiotemporal = 'rw', n.fold=5, tv=T,
                     x.formula.tv = xf3tv)


write.csv(rbind(mm.tv1,mm.tv2,mm.tv3 ), 'results/analysis_3d/pa_random_tv_performance.csv', row.names = F)

mm.tv1
summary(mm.tv1)
model.performance.s[9,]
mm.tv1

best.tv=sdmTMB(
  formula = pa ~ 0 + sal.fal+ tem.fal,
  time_varying = ~ 0 +  dep,
  #time_varying_type = 'rw',
  data = x.dat.m,
  mesh = mesh,
  family = binomial(link = "logit"),
  spatial='on',
  time = "year",
  spatiotemporal='rw'
)
saveRDS(best.tv, 'results/analysis_3d/pa_tv_model.RDS')

pred.grid=as.data.frame(basegrid)%>%
  dplyr::select(year, fyear,dep, sal.fal, tem.fal, FID)
pred.grid$Xkm=coords.grid$X/1000
pred.grid$Ykm=coords.grid$Y/1000
x.pred.grid=predict(best.tv, newdata=pred.grid)
pred.grid$pred.tv=best.tv$family$linkinv(x.pred.grid$est)
pred.grid=pred.grid[,c('year', 'FID','pred.tv')]

basegrid=left_join(basegrid, pred.grid, by=c('year','FID'))
basegrid$diff=basegrid$pred-basegrid$pred.tv

ggplot(data=basegrid)+
  geom_sf(aes(fill=diff), color=NA)+
  facet_wrap(~year)+
  scale_fill_viridis_c()


saveRDS(basegrid, 'results/analysis_3d/grid_pa.RDS')


summary(best.tv)


nd <- expand.grid(
  dep = seq(min(x.dat.m$dep) + 0.2,
                     max(x.dat.m$dep) - 0.2,
                     length.out = 50
  ),
  sal.fal=mean(x.dat.m$sal.fal),
  tem.fal=mean(x.dat.m$tem.fal),
  year = unique(x.dat.m$year)
)
#nd$fyear=as.factor(nd$year)


p <- predict(best.tv, newdata = nd, se_fit = TRUE, re_form = NA)
p$r.s=ifelse(p$year<2013,'regime1','regime2')
ggplot(p, aes(dep, exp(est),
             # ymin = exp(est - 1.96 * est_se),
             # ymax = exp(est + 1.96 * est_se),
              group = as.factor(year))) +
  geom_line(aes(colour = year, linetype=factor(r.s)), lwd = 1) +
  #geom_ribbon(aes(fill = year), alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Depth z.scored", y = "log recruit abundance density (n/km2)")

ggsave('results/analysis_3d/timevarying_dep.png', width = 20, height = 15, units='cm')

df.res=data.frame(tidy(best.tv, conf.int = TRUE))
p=ggplot(data=df.res[df.res$term!='(Intercept)',])+
  geom_vline(xintercept = 0)+
  geom_point(aes(y=term, x=estimate))+
  geom_errorbar(aes(y=term, xmin=conf.low, xmax=conf.high))
xpreds=simulate(best.tv, nsim = 999)
x.dat.m$pred.tv=apply(xpreds, 1, median)

caret::confusionMatrix(table(x.dat.m$pred, x.dat.m$pa))


x.dat.m$resid.tv=(apply(xpreds, 1, median))-x.dat.m$pa
ggplot(data=x.dat.m,aes(x=Xkm, y=Ykm))+
  geom_point(aes( pch=factor(pa), color=factor(resid.tv)), size=3)+
  facet_wrap(~year)

plot(x.dat.m$pred, x.dat.m$pa)
abline(a=0,b=1)

















