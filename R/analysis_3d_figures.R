rm(list=ls())
scriptPath <- rstudioapi::getSourceEditorContext()$path 
scriptDir <- dirname(scriptPath)
setwd(file.path(scriptDir, '..'))
library(tidyverse)
library(sf)
library(sdmTMB) # ???0.6.0???
library(spdep)
library(visreg)
library(viridis)


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

## Spatial predictions ####
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



### Figure 4
library(rnaturalearth)
library(ggspatial)
fig4.df=as.data.frame(basegrid)%>%
  #dplyr::filter(index<thr.out)%>%
  #dplyr::filter(index>thr)%>%
  dplyr::select(year,FID,index)%>%
  pivot_wider(names_from = year, values_from = index)%>%
  replace(is.na(.),0)%>%
  pivot_longer(-c(FID), names_to = 'year', values_to = 'index')%>%
  dplyr::mutate(regime=ifelse(year<=2012, 'Regime1','Regime2'))%>%
  dplyr::group_by(FID, regime)%>%
  dplyr::summarise(index=mean(index))%>%
  pivot_wider(names_from = regime, values_from = index)%>%
  dplyr::mutate(diff=log(Regime2+1)-log(Regime1+1))%>%
  dplyr::left_join(basegrid)%>%
  st_as_sf()
saveRDS(fig4.df, 'results/spatial_predictions.RDS')


spat.df=fig4.df%>%
  dplyr::group_by(x)%>%
  dplyr::summarise(Regime1=mean(Regime1), Regime2=mean(Regime2))%>%
  dplyr::mutate(Diff=Regime2-Regime1)
spat.smooth=st_intersection(spat.df, solemon.area)

countries=ne_countries(continent='europe', scale=10)%>%
  st_set_crs(4326)%>%
  st_transform(3003)%>%
  st_crop(st_buffer(spat.df, 20*1000))

ita=ne_countries(country='italy', scale=10)%>%
  st_set_crs(4326)%>%
  st_transform(3003)%>%
  st_crop(st_buffer(spat.df,15000))

closure3=st_buffer(ita, 4*1852)%>%
  st_difference(ita)%>%
  st_crop(solemon.area)

closure6=st_buffer(ita, 6*1852)%>%
  st_difference(ita)%>%
  st_crop(solemon.area)

range_vals <- range(log(c(spat.smooth$Regime1, spat.smooth$Regime2)+1), na.rm = TRUE)

set_theme(theme_bw())
p1=ggplot()+
  geom_sf(data=spat.smooth, aes(fill=log(Regime1+1)), color=NA)+
  geom_sf(data=closure3, fill=NA, color='red')+
  geom_sf(data=countries, fill='grey50')+
  labs(fill='Common sole juveniles (< 20 cm TL): log abundance (n/km2)')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank())+
  annotation_scale()+
  scale_fill_viridis(limits = range_vals)

p2=ggplot()+
  geom_sf(data=spat.smooth, aes(fill=log(Regime2+1)), color=NA)+
  geom_sf(data=closure6, fill=NA, color='red')+
  geom_sf(data=countries, fill='grey50')+
  labs(fill='log juveniles')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank())+
  annotation_scale()+
  scale_fill_viridis(limits = range_vals)

p3=ggplot()+
  geom_sf(data=spat.smooth, aes(fill=log(Diff+1)), color=NA)+
  geom_sf(data=closure6, fill=NA, color='red')+
  geom_sf(data=countries, fill='grey50')+
  labs(fill='log juveniles')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank())+
  annotation_scale()+
  scale_fill_viridis(limits = range_vals)

pcomb=ggpubr::ggarrange(p1,p2,p3, common.legend = T, nrow=1, labels=c('a)','b)','c)'))

ggsave(plot=pcomb, 'C:/github/solea_regime_shift/results/analysis_3d/plots/Fig4.jpeg', width = 20, height = 10, units='cm', dpi=500)

## annual plots

pl.yr=ggplot()+
  geom_sf(data=basegrid, aes(fill=log(index+1)), color=NA)+
  geom_sf(data=countries, fill='grey50')+
  labs(fill=' juvenile log abundance')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.position = 'bottom')+
  annotation_scale()+
  scale_fill_viridis(limits = range_vals)+
  facet_wrap(~year)

ggsave(plot=pl.yr, 'results/analysis_3d/plots/spatialyr_tv.jpeg', width = 20, height = 30, units='cm')



