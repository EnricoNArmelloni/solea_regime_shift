rm(list=ls())
scriptPath <- rstudioapi::getSourceEditorContext()$path 
scriptDir <- dirname(scriptPath)
setwd(file.path(scriptDir, '..'))
library(tidyverse)
library(sf)
library(sdmTMB) # ‘0.6.0’
library(spdep)
library(visreg)



# Table of model selection ####
mod.res=list.files('results', pattern='.csv')

# p/a
pa.mod=mod.res[grep('pa', mod.res)]
xres=NULL

pa.model= readRDS("results/analysis_3d/pa_model.RDS")
pa.tv.model= readRDS("results/analysis_3d/pa_tv_model.RDS")
pos.model= readRDS("results/analysis_3d/pos_model.RDS")
pos.tv.model= readRDS("results/analysis_3d/pos_tv_model.RDS")
xdat=read_csv("data/model_3d_input.csv")
m1.dat= read_csv("data/input_model_pa.csv")
m2.dat= read_csv("data/input_model_pos.csv")

## conditional effects ####
# PA: dep tv, sal, tem
pa.vec.dep <- expand.grid(
  dep = seq(min(m1.dat$dep),max(m1.dat$dep),length.out = 50),
  sal.fal=mean(m1.dat$sal.fal),
  tem.fal=mean(m1.dat$tem.fal),
  year = unique(m1.dat$year))
pa.vec.dep$depth=(pa.vec.dep$dep*sd(xdat$meandepth))+mean(xdat$meandepth)
p.dep <- predict(pa.tv.model, newdata = pa.vec.dep, se_fit = TRUE, re_form = NA)
p.dep$pred=pa.tv.model$family$linkinv(p.dep$est)
p.dep$r.s=ifelse(p.dep$year<2013,'regime1','regime2')


pa.vec.sal <- expand.grid(
  dep = mean(m1.dat$dep),
  sal.fal=seq(min(m1.dat$dep),max(m1.dat$dep),length.out = 50),
  tem.fal=mean(m1.dat$tem.fal),
  year = 2013)
pa.vec.sal$salinity=(pa.vec.sal$sal.fal*sd(xdat$sal_autumn))+mean(xdat$sal_autumn)
p.sal <- predict(pa.tv.model, newdata = pa.vec.sal, se_fit = TRUE, re_form = NA)
p.sal$pred=pa.tv.model$family$linkinv(p.sal$est)
p.sal$pred_se=pa.tv.model$family$linkinv(p.sal$est_se)


pa.vec.tem <- expand.grid(
  dep = mean(m1.dat$dep),
  sal.fal=mean(m1.dat$sal.fal),
  tem.fal=seq(min(m1.dat$dep),max(m1.dat$dep),length.out = 50),
  year = 2013)
pa.vec.tem$temperature=(pa.vec.tem$tem.fal*sd(xdat$bottomT_autumn))+mean(xdat$bottomT_autumn)
p.tem <- predict(pa.tv.model, newdata = pa.vec.tem, se_fit = TRUE, re_form = NA)
p.tem$pred=pa.tv.model$family$linkinv(p.tem$est)
p.tem$pred_se=pa.tv.model$family$linkinv(p.tem$est_se)


ppa1=ggplot(p.dep, aes(depth, pred,
              # ymin = exp(est - 1.96 * est_se),
              # ymax = exp(est + 1.96 * est_se),
              group = as.factor(year))) +
  geom_line(aes(colour = year, linetype=factor(r.s)), lwd = 1) +
  #geom_ribbon(aes(fill = year), alpha = 0.1) +
  scale_colour_viridis_c() +
  labs(linetype='Regime', colour='Year')+
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Depth", y = "Conditional effect on P/A")

ppa2=ggplot(p.sal, aes(salinity, pred,
                   ymin = (pred -  pred_se/2),
                   ymax = (pred +  pred_se/2),
                  group = as.factor(year))) +
  geom_line() +
  geom_ribbon( alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Salinity", y = "Conditional effect on P/A")

ppa3=ggplot(p.tem, aes(temperature, pred,
                  ymin = (pred -  pred_se/2),
                  ymax = (pred +  pred_se/2),
                  group = as.factor(year))) +
  geom_line() +
  geom_ribbon( alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Temperature", y = "Conditional effect on P/A")

ppa.comb=ggpubr::ggarrange(ppa1,ppa2,ppa3, nrow=1)

# pos: sal summer, depth, temfal tv
pos.vec.tem<- expand.grid(
  tem.fal = seq(min(m2.dat$tem.fal),max(m2.dat$tem.fal),length.out = 50),
  sal.sum=mean(m2.dat$sal.sum),
  ppr.sum=mean(m2.dat$ppr.sum),
  dep=mean(m2.dat$dep),
  year = unique(m2.dat$year))
pos.vec.tem$temperature=(pos.vec.tem$tem.fal*sd(xdat$bottomT_autumn))+mean(xdat$bottomT_autumn)
pos.tem <- predict(pos.tv.model, newdata = pos.vec.tem, se_fit = TRUE, re_form = NA)
pos.tem$r.s=ifelse(pos.tem$year<2013,'regime1','regime2')
pos.tem$pred=pos.tv.model$family$linkinv(pos.tem$est)

pos.vec.dep<- expand.grid(
  tem.fal = mean(m2.dat$tem.fal),
  sal.sum=mean(m2.dat$sal.sum),
  ppr.sum=mean(m2.dat$ppr.sum),
  dep=seq(min(m2.dat$dep),max(m2.dat$dep),length.out = 50),
  year = 2012)
pos.vec.dep$depth=(pos.vec.dep$dep*sd(xdat$meandepth))+mean(xdat$meandepth)
pos.dep <- predict(pos.tv.model, newdata = pos.vec.dep, se_fit = TRUE, re_form = NA)
pos.dep$r.s=ifelse(pos.dep$year<2013,'regime1','regime2')
pos.dep$pred=pos.tv.model$family$linkinv(pos.dep$est)
pos.dep$pred_se=pos.tv.model$family$linkinv(pos.dep$est_se)

pos.vec.sal<- expand.grid(
  tem.fal = mean(m2.dat$tem.fal),
  sal.sum=seq(min(m2.dat$sal.sum),max(m2.dat$sal.sum),length.out=50),
  ppr.sum=mean(m2.dat$ppr.sum),
  dep=mean(m2.dat$dep),
  year = 2012)
pos.vec.sal$salinity=(pos.vec.sal$sal.sum*sd(xdat$sal_summer))+mean(xdat$sal_summer)
pos.sal <- predict(pos.tv.model, newdata = pos.vec.sal, se_fit = TRUE, re_form = NA)
pos.sal$r.s=ifelse(pos.sal$year<2013,'regime1','regime2')
pos.sal$pred=pos.tv.model$family$linkinv(pos.sal$est)
pos.sal$pred_se=pos.tv.model$family$linkinv(pos.sal$est_se)

ppos1=ggplot(pos.tem, aes(temperature, pred,
                  # ymin = exp(est - 1.96 * est_se),
                  # ymax = exp(est + 1.96 * est_se),
                  group = as.factor(year))) +
  geom_line(aes(colour = year, linetype=factor(r.s)), lwd = 1) +
  #geom_ribbon(aes(fill = year), alpha = 0.1) +
  labs(linetype='Regime', colour='Year')+
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Temperature", y = "Conditional effect on Density")

ppos2=ggplot(pos.dep, aes(depth, pred,
                          ymin = (pred -  pred_se/2),
                          ymax = (pred +  pred_se/2),
                    group = as.factor(year))) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Depth", y = "Conditional effect on Density")

ppos3=ggplot(pos.sal, aes(salinity, pred,
                          ymin = (pred -  pred_se/2),
                          ymax = (pred +  pred_se/2),
                    group = as.factor(year))) +
  geom_line() +
  geom_ribbon( alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  labs(x = "Salinity", y = "Conditional effect on Density")


pcomb=ggpubr::ggarrange(ppa1,ppa2,ppa3,ppos2,ppos3, ppos1, nrow=2, ncol=3 ,common.legend = T)


#ptv3=ggpubr::ggarrange(ptv1, ptv2, common.legend = T)
ggsave(plot=pcomb,'results/analysis_3d/plots/conditional_effects_tv.jpeg', width = 15, height = 10, units='cm')


## residuals ####
df.res=data.frame(tidy(pa.model, conf.int = TRUE))
p=ggplot(data=df.res[df.res$term!='(Intercept)',])+
  geom_vline(xintercept = 0)+
  geom_point(aes(y=term, x=estimate))+
  geom_errorbar(aes(y=term, xmin=conf.low, xmax=conf.high))

# residuals plots ##
# qq plot
set.seed(123)
rq_res <- residuals(pa.model, type = "mle-mvn")
rq_res <- rq_res[is.finite(rq_res)] # some Inf

xpreds=simulate(pa.model)
m1.dat$pred=apply(xpreds, 1, median)

m1.dat$residuals=rq_res
m1.dat$resid=(apply(xpreds, 1, median))-m1.dat$pa


jpeg('results/analysis_3d/plots/res_pa.jpeg', width = 15, height = 10, units='cm', res=200)
par(mfrow=c(1,2))
qqnorm(rq_res);abline(0, 1)
hist(rq_res)
dev.off()

ggplot(data=m1.dat,aes(x=Xkm, y=Ykm))+
  geom_point(aes( pch=factor(pa), color=factor(resid)), size=3)+
  facet_wrap(~year)

pr1=ggplot(data=m1.dat, aes(x=dep, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr2=ggplot(data=m1.dat, aes(x=sal.fal, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr3=ggplot(data=m1.dat, aes(x=tem.fal, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr4=ggplot(data=m1.dat, aes(x=year, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)

p12=ggpubr::ggarrange(pr1, pr2)
p34=ggpubr::ggarrange(pr3, pr4)
p5=ggpubr::ggarrange(p12, p34, ncol=1)
ggsave(plot=p5, 'results/analysis_3d/plots/residuals_pa.png', width=15, height=15, units='cm')

# pa tv
set.seed(123)
rq_res <- residuals(pa.tv.model, type = "mle-mvn")
rq_res <- rq_res[is.finite(rq_res)] # some Inf

xpreds=simulate(pa.tv.model)
m1.dat$pred=apply(xpreds, 1, median)

m1.dat$residuals=rq_res
m1.dat$resid=(apply(xpreds, 1, median))-m1.dat$pa


jpeg('results/analysis_3d/plots/res_pa_tv.jpeg', width = 15, height = 10, units='cm', res=200)
par(mfrow=c(1,2))
qqnorm(rq_res);abline(0, 1)
hist(rq_res)
dev.off()

ggplot(data=m1.dat,aes(x=Xkm, y=Ykm))+
  geom_point(aes( pch=factor(pa), color=factor(resid)), size=3)+
  facet_wrap(~year)

pr1=ggplot(data=m1.dat, aes(x=dep, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr2=ggplot(data=m1.dat, aes(x=sal.fal, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr3=ggplot(data=m1.dat, aes(x=tem.fal, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr4=ggplot(data=m1.dat, aes(x=year, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)

p12=ggpubr::ggarrange(pr1, pr2)
p34=ggpubr::ggarrange(pr3, pr4)
p5=ggpubr::ggarrange(p12, p34, ncol=1)
ggsave(plot=p5, 'results/analysis_3d/plots/residuals_pa_tv.png', width=15, height=15, units='cm')

# pos
set.seed(123)
rq_res <- residuals(pos.model, type = "mle-mvn")
rq_res <- rq_res[is.finite(rq_res)] # some Inf

xpreds=simulate(pos.model)
m2.dat$pred=apply(xpreds, 1, median)

m2.dat$residuals=rq_res
m2.dat$resid=(apply(xpreds, 1, median))-m2.dat$logrec



jpeg('results/analysis_3d/plots/res_pos.jpeg', width = 15, height = 10, units='cm', res=200)
par(mfrow=c(1,2))
qqnorm(rq_res);abline(0, 1)
hist(rq_res)
dev.off()

ggplot(data=m2.dat,aes(x=Xkm, y=Ykm))+
  geom_point(aes( color=(resid)), size=3)+
  facet_wrap(~year)

pr1=ggplot(data=m2.dat, aes(x=dep, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr2=ggplot(data=m2.dat, aes(x=sal.sum, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr3=ggplot(data=m2.dat, aes(x=ppr.sum, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr4=ggplot(data=m2.dat, aes(x=tem.fal, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr5=ggplot(data=m2.dat, aes(x=year, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)

p12=ggpubr::ggarrange(pr1, pr2,nrow=1)
p34=ggpubr::ggarrange(pr4, pr5)
p5=ggpubr::ggarrange(p12, p34, ncol=1)
ggsave(plot=p5, 'results/analysis_3d/plots/residuals_pos.png', width=15, height=15, units='cm')

# pos tv
set.seed(123)
rq_res <- residuals(pos.tv.model, type = "mle-mvn")
rq_res <- rq_res[is.finite(rq_res)] # some Inf

xpreds=simulate(pos.tv.model)
m2.dat$pred=apply(xpreds, 1, median)

m2.dat$residuals=rq_res
m2.dat$resid=(apply(xpreds, 1, median))-m2.dat$logrec

qqnorm(log(m2.dat$resid))


jpeg('results/analysis_3d/plots/res_pos_tv.jpeg', width = 15, height = 10, units='cm', res=200)
par(mfrow=c(1,2))
qqnorm(rq_res);abline(0, 1)
hist(rq_res)
dev.off()

ggplot(data=m2.dat,aes(x=Xkm, y=Ykm))+
  geom_point(aes( color=(resid)), size=3)+
  facet_wrap(~year)

pr1=ggplot(data=m2.dat, aes(x=dep, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr2=ggplot(data=m2.dat, aes(x=sal.sum, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr3=ggplot(data=m2.dat, aes(x=ppr.sum, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr4=ggplot(data=m2.dat, aes(x=tem.fal, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)
pr5=ggplot(data=m2.dat, aes(x=year, y=residuals))+
  geom_point()+
  geom_hline(yintercept = 0)

p12=ggpubr::ggarrange(pr1, pr2, nrow=1)
p34=ggpubr::ggarrange(pr4, pr5)
p5=ggpubr::ggarrange(p12, p34, ncol=1)
ggsave(plot=p5, 'results/analysis_3d/plots/residuals_pos_tv.png', width=15, height=15, units='cm')


## combined predictions residuals
m1.dat$fyear=as.factor(m1.dat$fyear)
m1.dat$rec=exp(m1.dat$logrec)-1
m1.dat$sqrec=sqrt(m1.dat$rec)

pred.pa=predict(pa.model, newdata=m1.dat)
pred.pa=pa.model$family$linkinv(pred.pa$est)
pred.pos=predict(pos.model, newdata=m1.dat)
pred.pos=pos.model$family$linkinv(pred.pos$est)
pred.comb=pred.pa*pred.pos

m1.dat$prediction.comb.base=pred.comb

pred.pa=predict(pa.tv.model, newdata=m1.dat)
pred.pa=pa.tv.model$family$linkinv(pred.pa$est)
pred.pos=predict(pos.tv.model, newdata=m1.dat)
pred.pos=pos.tv.model$family$linkinv(pred.pos$est)
pred.comb=pred.pa*pred.pos

m1.dat$prediction.comb.tv=pred.comb
m1.dat$res.base=(m1.dat$sqrec - m1.dat$prediction.comb.base)
m1.dat$res.tv=(m1.dat$sqrec - m1.dat$prediction.comb.tv)


jpeg('results/analysis_3d/plots/pred_comb.jpeg', width = 25, height = 10, units='cm', res=200)
par(mfrow=c(1,2))
plot(m1.dat$sqrec, m1.dat$prediction.comb.base, main='base')
abline(a=0,b=1)
plot(m1.dat$sqrec, m1.dat$prediction.comb.tv, main='tv')
abline(a=0,b=1)
dev.off()

ggplot(data=m2.dat,aes(x=Xkm, y=Ykm))+
  geom_point(aes( color=(resid)), size=3)+
  facet_wrap(~year)

pr1=ggplot(data=m1.dat, aes(x=dep, y=res.base))+
  geom_point()+
  geom_hline(yintercept = 0)
pr2=ggplot(data=m1.dat, aes(x=sal.sum, y=res.base))+
  geom_point()+
  geom_hline(yintercept = 0)
pr3=ggplot(data=m1.dat, aes(x=ppr.sum, y=res.base))+
  geom_point()+
  geom_hline(yintercept = 0)
pr4=ggplot(data=m1.dat, aes(x=tem.fal, y=res.base))+
  geom_point()+
  geom_hline(yintercept = 0)
pr5=ggplot(data=m1.dat, aes(x=year, y=res.base))+
  geom_point()+
  geom_hline(yintercept = 0)

p12=ggpubr::ggarrange(pr1, pr2, nrow=1)
p34=ggpubr::ggarrange(pr4, pr5)
p5=ggpubr::ggarrange(p12, p34, ncol=1)
ggsave(plot=p5, 'results/analysis_3d/plots/residuals_base.png', width=15, height=15, units='cm')


pr1=ggplot(data=m1.dat, aes(x=dep, y=res.tv))+
  geom_point()+
  geom_hline(yintercept = 0)
pr2=ggplot(data=m1.dat, aes(x=sal.sum, y=res.tv))+
  geom_point()+
  geom_hline(yintercept = 0)
pr3=ggplot(data=m1.dat, aes(x=ppr.sum, y=res.tv))+
  geom_point()+
  geom_hline(yintercept = 0)
pr4=ggplot(data=m1.dat, aes(x=tem.fal, y=res.tv))+
  geom_point()+
  geom_hline(yintercept = 0)
pr5=ggplot(data=m1.dat, aes(x=year, y=res.tv))+
  geom_point()+
  geom_hline(yintercept = 0)

p12=ggpubr::ggarrange(pr1, pr2, nrow=1)
p34=ggpubr::ggarrange(pr4, pr5)
p5=ggpubr::ggarrange(p12, p34, ncol=1)
ggsave(plot=p5, 'results/analysis_3d/plots/residuals_tv.png', width=15, height=15, units='cm')



#s_pa <- simulate(best.model.pa, nsim = 500, type = "mle-mvn")
#r_pa=dharma_residuals(s_pa, best.model.pa, return_DHARMa = T)
#DHARMa::testResiduals(r_pa)
#DHARMa::testSpatialAutocorrelation(r_pa, x = x.dat.m$Xkm, y = x.dat.m$Ykm)
#library(DHARMa)


caret::confusionMatrix(table(x.dat.m$pred, x.dat.m$pa))

# spatial preds ####
#solemon.area=read_sf('C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/Data/shapefiles/Solemon_strata_ITA_SVN_depth')%>%
solemon.area=read_sf('../other_data/Solemon_strata_ITA_SVN_depth')%>%
  st_union()%>%
  st_set_crs(4326)%>%
  st_transform(3003)

basegrid <- readRDS("data/basegrid.RDS")
basegrid$depth=-basegrid$depth
basegrid$dep=(basegrid$depth-mean(xdat$meandepth))/sd(xdat$meandepth)
basegrid$sal.fal=(basegrid$sal_autumn-mean(xdat$sal_autumn))/sd(xdat$sal_autumn)
basegrid$tem.fal=(basegrid$bottomT_autumn-mean(xdat$bottomT_autumn))/sd(xdat$bottomT_autumn)
basegrid$fyear=as.factor(basegrid$year)
basegrid$ppr.sum=(basegrid$pp_summer-mean(xdat$pp_summer))/sd(xdat$pp_summer)
basegrid$sal.sum=(basegrid$sal_summer-mean(xdat$sal_summer))/sd(xdat$sal_summer)

#med=read_sf("C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/BigData/contours/Med_Poly")%>%
med=read_sf("../other_data/Med_Poly")%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)%>%
  st_crop(.,st_buffer(basegrid))

coords.grid=as.data.frame(st_coordinates(st_centroid(basegrid)))

pred.grid=as.data.frame(basegrid)%>%
  dplyr::select(year, fyear,dep, sal.fal, tem.fal, FID)
pred.grid$Xkm=coords.grid$X/1000
pred.grid$Ykm=coords.grid$Y/1000

x.pred.grid=predict(pa.model, newdata=pred.grid)
pred.grid$pred=pa.model$family$linkinv(x.pred.grid$est)
pred.grid$rf.pa=x.pred.grid$omega_s
pred.grid=pred.grid[,c('year', 'FID','pred', 'rf.pa')]

basegrid=left_join(basegrid, pred.grid, by=c('year','FID'))

pred.grid=as.data.frame(basegrid)%>%
  dplyr::select(year, fyear,dep, sal.fal,sal.sum, ppr.sum, tem.fal,FID)
pred.grid$Xkm=coords.grid$X/1000
pred.grid$Ykm=coords.grid$Y/1000
x.pred.grid=predict(pos.model, newdata=pred.grid)
pred.grid$pred.pos=pos.model$family$linkinv(x.pred.grid$est)
pred.grid=pred.grid[,c('year', 'FID','pred.pos')]

basegrid=left_join(basegrid, pred.grid, by=c('year','FID'))
basegrid$index=(basegrid$pred.pos)*(basegrid$pred^2)

#basegrid$index.tv=(basegrid$pred.pos)*basegrid$pred.tv

#saveRDS(basegrid, 'results/grid_pos.RDS')
thr.out=as.numeric(quantile(basegrid$index, probs=0.999))

pl.yr=ggplot(data=basegrid)+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(aes(fill=log(index+1)), color=NA)+
  facet_wrap(~year)+
  scale_fill_viridis_c()+
  geom_sf(data=med)

ggsave(plot=pl.yr, 'results/analysis_3d/plots/spatialyr.jpeg', width = 20, height = 30, units='cm')

thr=as.numeric(quantile(basegrid$index, probs=0.95))

test=basegrid%>%
  #dplyr::filter(index<thr.out)%>%
  dplyr::filter(index>thr)%>%
  dplyr::mutate(regime=ifelse(year<2011, 'Regime 1','Regime 2 (>2010)'))%>%
  dplyr::group_by(FID, regime)%>%
  dplyr::summarise(index=mean(index))%>%
  dplyr::group_by(FID)%>%
  dplyr::mutate(diff=index-mean(index))
  
p.persistence=ggplot(data=test)+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(aes(fill=log(index+1)), color=NA)+
  facet_wrap(~regime)+
  scale_fill_viridis_c()+
  geom_sf(data=med)+
  labs(fill='log recruits (n/km2)')+
  theme(legend.position = 'bottom')

ggsave(plot=p.persistence, 'plots/persistence.jpeg', width = 10, height = 10, units='cm')


# spatial preds tv ####
solemon.area=read_sf('../other_data/Solemon_strata_ITA_SVN_depth')%>%
  st_union()%>%
  st_set_crs(4326)%>%
  st_transform(3003)

basegrid <- readRDS("data/basegrid.RDS")
basegrid$depth=-basegrid$depth
basegrid$dep=(basegrid$depth-mean(xdat$meandepth))/sd(xdat$meandepth)
basegrid$sal.fal=(basegrid$sal_autumn-mean(xdat$sal_autumn))/sd(xdat$sal_autumn)
basegrid$tem.fal=(basegrid$bottomT_autumn-mean(xdat$bottomT_autumn))/sd(xdat$bottomT_autumn)
basegrid$fyear=as.factor(basegrid$year)
basegrid$ppr.sum=(basegrid$pp_summer-mean(xdat$pp_summer))/sd(xdat$pp_summer)
basegrid$sal.sum=(basegrid$sal_summer-mean(xdat$sal_summer))/sd(xdat$sal_summer)

coords.grid=as.data.frame(st_coordinates(st_centroid(basegrid)))

pred.grid=as.data.frame(basegrid)%>%
  dplyr::select(year, fyear,dep, sal.fal, tem.fal, FID)
pred.grid$Xkm=coords.grid$X/1000
pred.grid$Ykm=coords.grid$Y/1000

x.pred.grid=predict(pa.tv.model, newdata=pred.grid)
pred.grid$pred=pa.tv.model$family$linkinv(x.pred.grid$est)
pred.grid$rf.pa=x.pred.grid$omega_s
pred.grid=pred.grid[,c('year', 'FID','pred', 'rf.pa')]

basegrid=left_join(basegrid, pred.grid, by=c('year','FID'))

pred.grid=as.data.frame(basegrid)%>%
  dplyr::select(year, fyear,dep, sal.fal,sal.sum, ppr.sum, tem.fal,FID)
pred.grid$Xkm=coords.grid$X/1000
pred.grid$Ykm=coords.grid$Y/1000
x.pred.grid=predict(pos.tv.model, newdata=pred.grid)
pred.grid$pred.pos=pos.tv.model$family$linkinv(x.pred.grid$est)
pred.grid=pred.grid[,c('year', 'FID','pred.pos')]

basegrid=left_join(basegrid, pred.grid, by=c('year','FID'))
basegrid$index=(basegrid$pred.pos)*(basegrid$pred)

#basegrid$index.tv=(basegrid$pred.pos)*basegrid$pred.tv

#saveRDS(basegrid, 'results/grid_pos.RDS')
thr.out=as.numeric(quantile(basegrid$index, probs=0.999))

base.stats=basegrid%>%dplyr::group_by(year)%>%
  dplyr::summarise(pa=mean(pred),pa.var=sd(pred), pos=mean(pred.pos), pos.var=sd(pred.pos))

ppasp=ggplot(data=base.stats)+
  geom_line(aes(x=year, y=pa))+
  geom_ribbon(aes(x=year,ymin=pa-pa.var/2, ymax=pa+pa.var/2), alpha=0.1)+
  geom_vline(xintercept=2012)+
  xlab('Year')+
  ylab('Average presence over study area')

ppossp=ggplot(data=base.stats)+
  geom_line(aes(x=year, y=pos))+
  geom_ribbon(aes(x=year,ymin=pos-pos.var/2, ymax=pos+pos.var/2), alpha=0.1)+
  geom_vline(xintercept=2012)+
  xlab('Year')+
  ylab('Average density over study area')

pcomb=ggpubr::ggarrange(ppa1,ppa2,ppa3,ppasp,ppos2,ppos3, ppos1,ppossp, nrow=2, ncol=4 ,common.legend = T)
ggsave(plot=pcomb,'results/analysis_3d/plots/conditional_effects_tv.jpeg', width = 25, height = 15, units='cm')



pl.yr=ggplot(data=basegrid)+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(aes(fill=index), color=NA)+
  facet_wrap(~year)+
  scale_fill_viridis_c()+
  geom_sf(data=med)

ggsave(plot=pl.yr, 'results/analysis_3d/plots/spatialyr_tv.jpeg', width = 20, height = 30, units='cm')

thr=as.numeric(quantile(basegrid$index, probs=0.95))

test=as.data.frame(basegrid)%>%
  #dplyr::filter(index<thr.out)%>%
  #dplyr::filter(index>thr)%>%
  dplyr::select(year,FID,index)%>%
  pivot_wider(names_from = year, values_from = index)%>%
  replace(is.na(.),0)%>%
  pivot_longer(-c(FID), names_to = 'year', values_to = 'index')%>%
  dplyr::mutate(regime=ifelse(year<2011, 'Regime1','Regime2'))%>%
  dplyr::group_by(FID, regime)%>%
  dplyr::summarise(index=mean(index))%>%
  pivot_wider(names_from = regime, values_from = index)%>%
  dplyr::mutate(diff=log(Regime2+1)-log(Regime1+1))%>%
  dplyr::left_join(basegrid)%>%
  st_as_sf()

saveRDS(test, 'results/spatial_predictions.RDS')

pm1=ggplot(data=test[test$Regime1>0,])+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(aes(fill=log(Regime1+1)), color=NA)+
  scale_fill_viridis_c()+
  geom_sf(data=med)+
  ggtitle('Regime 1')+
  labs(fill='log recruits (n/km2)')+
  theme(legend.position = 'bottom')

pm2=ggplot(data=test[test$Regime2>0,])+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(aes(fill=log(Regime2+1)), color=NA)+
  scale_fill_viridis_c()+
  geom_sf(data=med)+
  labs(fill='log recruits (n/km2)')+
  ggtitle('Regime 2')+
  theme(legend.position = 'bottom')

pm3=ggplot()+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(data=test, color=NA,
          aes(fill=(diff)))+
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0)+
  labs(fill='R2 - R1')+
  geom_sf(data=med)+
  ggtitle('Regime 2 - Regime 1')+
  theme(legend.position = 'bottom')

pr4=ggpubr::ggarrange(pm1,pm2, pm3, nrow=1)


ggsave(plot=pr4, 'results/analysis_3d/plots/persistence_tv.jpeg', width = 30, height = 15, units='cm')



t1=test[test$Regime1>0,]%>%
  dplyr::group_by(x)%>%
  dplyr::summarise(Regime1=mean(Regime1))
t2=test[test$Regime2>0,]%>%
  dplyr::group_by(x)%>%
  dplyr::summarise(Regime1=mean(Regime2))

pt1=ggplot(data=t1)+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(aes(fill=log(Regime1+1)), color=NA)+
  scale_fill_viridis_c()+
  geom_sf(data=med)+
  geom_sf(data=closure, fill=NA, color='red')+
  ggtitle('Regime 1')+
  labs(fill='log recruits (n/km2)')+
  theme(legend.position = 'bottom')

pt2=ggplot(data=t2)+
  geom_sf(data=solemon.area, fill='white')+
  geom_sf(aes(fill=log(Regime1+1)), color=NA)+
  scale_fill_viridis_c()+
  geom_sf(data=med)+
  geom_sf(data=closure, fill=NA, color='red')+
  ggtitle('Regime 1')+
  labs(fill='log recruits (n/km2)')+
  theme(legend.position = 'bottom')


ggpubr::ggarrange(pt1, pt2)

### hot spot
library(spdep)
hslist=list()

for(j in 1:length(unique(basegrid$year))){
  hsdat=basegrid[basegrid$year==unique(basegrid$year)[j], ]%>%st_as_sf()
  #hsdat=st_as_sf(hsdat, coords=c('Xkm','Ykm'))
  #hslist[[i]]=hsdat[hsdat$finalpred>quantile(hsdat$finalpred, probs = seq(0, 1, 0.1))[10], ]
  hsdat=na.omit(hsdat)
  nb <- dnearneigh(st_centroid(hsdat), d1=(1.852*2.5)*1000, d2=(1.852*5)*1000)
  nb_lw <- nb2listw(nb, zero.policy=T)
  local_g <- localG(hsdat$index, nb_lw)
  hsdat$local_g=as.numeric(local_g)
  hsdat=na.omit(hsdat)
  #plot(hsdat$x)
  #hsdat=hsdat[hsdat$local_g> as.numeric(quantile(hsdat$local_g, probs = seq(0, 1, 0.05)))[19], ]
  hslist[[j]]=hsdat
}
hslist=plyr::ldply(hslist)
hslist=hslist[hslist$local_g> as.numeric(quantile(hslist$local_g, probs = 0.9)), ]


r1=2007:2010
r2=2011:2019

hsresult=data.frame(hslist)%>%dplyr::select(FID, year)%>%
  dplyr::filter(year%in%r2)%>%
  dplyr::group_by(FID)%>%
  tally(name='hs')%>%
  dplyr::mutate(hs_class=(hs/length(r2))*100,
                hs_cat_p2=as.numeric(cut(hs_class, breaks=c(-Inf, 20,40,60,80,Inf), labels=1:5)))%>%
  dplyr::select(FID,hs_cat_p2 )
hsresultp2=na.omit(hsresult)

hsresult=data.frame(hslist)%>%dplyr::select(FID, year)%>%
  dplyr::filter(year%in%r1)%>%
  dplyr::group_by(FID)%>%
  tally(name='hs')%>%
  dplyr::mutate(hs_class=(hs/length(r1))*100,
                hs_cat_p1=as.numeric(cut(hs_class, breaks=c(-Inf, 20,40,60,80,Inf), labels=1:5)))%>%
  dplyr::select(FID,hs_cat_p1 )
hsresultp1=na.omit(hsresult)


hsresultcomb=left_join(basegrid,hsresultp2)%>%
  full_join(hsresultp1)%>%
  replace(is.na(.),0)%>%st_as_sf()
hsresultcomb$diff=hsresultcomb$hs_cat_p2-hsresultcomb$hs_cat_p1

library(rnaturalearth)
library(sf)
library(ggspatial)
xcountry=ne_countries(country = c("italy", 'croatia', 'slovenia', 'bosnia'), scale = "medium")%>%
  st_as_sf%>%
  st_set_crs(4326)%>%
  st_transform(3003)

med=read_sf("../other_data/Med_Poly")%>%
  st_set_crs(4326)%>%
  st_transform(., 3003)%>%
  st_crop(.,st_buffer(basegrid, 20000))


pl.map=ggplot()+
  geom_sf(data=med, fill='gray50')+
  geom_sf(data=solemon.area, fill='white', alpha=0.4)+
  annotation_scale() +
  labs(fill='Hot Spot category')+
  #theme_bw()+
  theme(panel.background = element_rect(fill='white'))+
  scale_fill_viridis_d()+
  #geom_sf(data=med)+
  theme(legend.position = 'bottom');pl.map

mycols=viridis(11)
mycols <- c(
  "#000000", # black
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "white", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # grey
  "#66C2A5", # teal-ish
  "#FC8D62"  # coral
)
names(mycols)=as.character(-5:5)

pr2=pl.map+geom_sf(data=hsresultcomb, color=NA,
                   aes(fill=factor(hs_cat_p2)))+
  scale_fill_manual(values = mycols);pr2



pr1=pl.map+geom_sf(data=hsresultcomb, color=NA,
                   aes(fill=factor(hs_cat_p1)))+
  scale_fill_manual(values = mycols);pr1

hsresultcomb$diffbin=ifelse(hsresultcomb$diff>0,1,
                            ifelse(hsresultcomb$diff<0,-1,0))
pr3=pl.map+
  #geom_sf(data=solemon.area, fill='black', alpha=0.4)+
  geom_sf(data=hsresultcomb, color=NA,
          aes(fill=as.factor(diff)))+
  scale_fill_manual(values = mycols)+
  #scale_fill_gradient2(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF", midpoint = 0)+
  labs(fill='Hotspot category changes');pr3

pr4=ggpubr::ggarrange(pr3,pr1,pr2,  nrow=1, common.legend = T);pr4

ggsave(plot=pr4, 'results/analysis_3d/plots/HotSpot_tv.jpeg', width = 30, height = 15, units='cm', dpi=300)


ggpubr::ggarrange(pr1,pr2, pr3, nrow=1, labels=c('a) R1: 2007-2012', 'b) R2: 2019-2019', 'c) R2 - R1'))


