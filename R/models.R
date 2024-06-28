setwd("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift")
rm(list=ls())
library(tidyverse)
library(sf)
library(sdmTMB)
library(spdep)

# load data
xdat=read.csv('data/model_input.csv')
land=read_sf("C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/contours/Med_Poly")%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)

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
head(xdat)
hist(log(xdat$Abun.rec+1))
hist(log(xdat$pp_spring))
hist(log(xdat$bottomT_autumn))

hist(log(xdat$hsi))
hist((xdat$grainsize))
hist(xdat$temp)
hist((xdat$sal_autumn))
hist(xdat$po_autumn)
hist(log(xdat$fhours_2019_OTB_autumn))
hist(log(xdat$fhours_2019_TBB_autumn))

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

xdat$pp.trans=transf.fun(xdat$pp_spring, 'z.score')
xdat$T.trans=transf.fun(xdat$bottomT_autumn, 'z.score')
xdat$dept.tranf=transf.fun(xdat$meandepth, 'z.score')

xdat.model$id=seq(1:nrow(xdat.model))

# look for spatial autocorrelation
xdatsf2=xdat%>%
  dplyr::group_by(xlon,xlat,year)%>%
  dplyr::summarise(rec=mean(Abun.rec))%>%
  st_as_sf(coords=c('xlon','xlat'))%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)
xdatsf2$logrec=log(xdatsf2$rec+1)

land.2=st_crop(land,st_buffer(xdatsf2,50000))
p.full=ggplot()+
  geom_sf(data=land.2, fill='darkgrey')+
  geom_sf(data=st_buffer(xdatsf2, 15000), aes(fill=logrec))+
  scale_fill_gradient2()+
  labs(fill='log recruitment density')+
  theme(legend.position = 'bottom')+
  facet_wrap(~year);p.full
ggsave('plots/recruit_distribution.png', width = 30, height = 50, units='cm')

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

## test some models
mod.res=function(x.formula){
  x.model <- sdmTMB(
    formula = x.formula,
    data = xdat,
    mesh = mesh,
    family = tweedie(link = "log"),
    spatial='on',
    time = "year",
    spatiotemporal='RW'
  ) # conv
  print(AIC(x.model))
  print(tidy(x.model, conf.int=T, 'ran_pars'))
  sanity(x.model)
  df.res=data.frame(tidy(x.model, conf.int = TRUE))
  ggplot(data=df.res[df.res$term!='(Intercept)',])+
    geom_vline(xintercept = 0)+
    geom_point(aes(y=term, x=estimate))+
    geom_errorbar(aes(y=term, xmin=conf.low, xmax=conf.high))
  
  xpreds=simulate(best.model, nsim = 999)
  xdat$pred=apply(xpreds, 1, median)
  xdat$resid=(apply(xpreds, 1, median))-xdat$logrec
  plot(xdat$pred, xdat$logrec)
  abline(a=0,b=1)
}

# endogenous variables
names(xdat.model)
xformula0=logrec ~ dept.tranf 
mod.res(xformula0) 
xformula0=logrec ~ pp.trans 
mod.res(xformula0) 
xformula0=logrec ~ T.trans
mod.res(xformula0) 
xformula0=logrec ~ dept.tranf + pp.trans  
mod.res(xformula0) 
xformula0=logrec ~ dept.tranf + T.trans
mod.res(xformula0) 
xformula0=logrec ~ dept.tranf + T.trans + pp.trans
mod.res(xformula0) 


# test time varying
names(xdat)
x.model <- sdmTMB(
  formula = logrec ~  0 ,
  data = xdat,
  time_varying = ~ 0 +  dept.tranf,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial='on',
  time = "year",
  spatiotemporal='ar1',
  silent=F
) 
x.model
AIC(x.model)
sanity(x.model)

nd <- expand.grid(
  dept.tranf = seq(min(xdat$dept.tranf) + 0.2,
                     max(xdat$dept.tranf) - 0.2,
                     length.out = 50
  ),
  year = unique(xdat$year) # all years
)

p <- predict(x.model, newdata = nd, se_fit = TRUE, re_form = NA)
p$r.s=ifelse(p$year<2013,'regime1','regime2')
ggplot(p, aes(dept.tranf, exp(est),
              ymin = exp(est - 1.96 * est_se),
              ymax = exp(est + 1.96 * est_se),
              group = as.factor(year)
)) +
  geom_line(aes(colour = year, linetype=factor(r.s)), lwd = 1) +
  #geom_ribbon(aes(fill = year), alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Depth z.scored", y = "log recruit abundance density (n/km2)")
ggsave('plots/timevarying_temp.png', width = 20, height = 15, units='cm')


