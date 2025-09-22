#########R flow to analyze Change Point and regime shifts in Solea Solea from Sguotti et al., ############
load("solea_tot.RData")#these are the data already prepared to perform the analyses and aggregated from the SoleMON data. 
#@AbunRecruits --> abundance recruits (fish smaller than 20cm) (thousands)
#@bottom_t --> Sea Bottom Temperature (C)
#@effort --> as calculated, then the values were standardized between 0 and 1. 
#@PP --> Primary Production
#@runoff_ms --> runoff river Po
#@AbunAdults --> Abundance adult fish

#Create a column with the total abundance of Sole to analyse it as a whole. 
solea_tot$Abun_tot <- solea_tot$AbunRecruits+solea_tot$AbunAdults
solea_tot$Abun_tot <- area_index_all$AbunIndex[1:13]


#####open packages that we will use for the analyses 
library(ggplot2)#plot
library(dplyr)
library(reshape2)
library(tidyr)
library(bcp)#change point Bayesian Change Point Analysis 
library(changepoint)#change point 

#my theme for plotting
cami_theme <- function() {
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(colour = "black", family = "Helvetica", size=8),
    plot.title = element_text(colour = "black", face = "bold", size = 10, hjust = 0.5, family = "Helvetica"),
    axis.title = element_text(colour = "black", face = "bold", size = 8, family = "Helvetica"),
    #panel.grid.major.x = element_line(colour = "dodgerblue3"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(family = "Helvetica", colour = "white", size=12),
    strip.background = element_rect(fill = "black"),
    axis.ticks = element_line(colour = "black")
  )
}


####1) plot the data #######
#plot recruitment over time 
p <- ggplot()
p <- p+ geom_point(data = solea_tot, aes(x =year, y =AbunRecruits ),color = "seagreen3",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =AbunRecruits ),color = "seagreen3")
p <- p+ ylab("Recruits (thousands)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("Recruits_sole.png", p, width=10, height=10, dpi=300,units="cm" )

#plot adults over time 
p <- ggplot()
p <- p+ geom_point(data = solea_tot, aes(x =year, y =AbunAdults ),color = "#41b6c4",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =AbunAdults ),color = "#41b6c4")
p <- p+ ylab("Adults (thousands)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("Adults_sole.png", p, width=10, height=10, dpi=300,units="cm" )

#plot total abundance over time 
p <- ggplot()
p <- p+ geom_point(data = solea_tot, aes(x =year, y =Abun_tot ),color = "#225ea8",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =Abun_tot ),color = "#225ea8")
p <- p+ ylab("Sole (number/sqkm)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("Tot_sole.png", p, width=10, height=10, dpi=300,units="cm" )

#time series temperature
p <- ggplot()
p <- p+ geom_point(data = solea_tot, aes(x =year, y =bottom_t ),color = "black",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =bottom_t),color = "black")
p <- p+ ylab("Sea Bottom Temperature (C)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("SBT.png", p, width=8, height=8, dpi=300,units="cm" )

#time series effort
p <- ggplot()
p <- p+ geom_point(data = solea_tot, aes(x =year, y =effort ),color = "black",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =effort),color = "black")
p <- p+ ylab("Fishing Capacity")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("Fishing.png", p, width=8, height=8, dpi=300,units="cm" )

#time series PP
p <- ggplot()
p <- p+ geom_point(data = solea_tot, aes(x =year, y =PP ),color = "black",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =PP),color = "black")
p <- p+ ylab("Primary Production (C per unit Volume)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("PP.png", p, width=8, height=8, dpi=300,units="cm" )

#time series runoff Po river 
p <- ggplot()
p <- p+ geom_point(data = solea_tot, aes(x =year, y =runoff_ms ),color = "black",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =runoff_ms),color = "black")
p <- p+ ylab("Runoff River Po (m³/s)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("Po.png", p, width=8, height=8, dpi=300,units="cm" )

################Change Point Analysis############

#change point 


cpt.mean(solea_tot$AbunRecruits, Q=1,method="BinSeg")#location4


cpt.mean(solea_tot$AbunAdults,Q=1,method="BinSeg")#location 7


cpt.mean(solea_tot$Abun_tot,Q=1,method="BinSeg")#location 6


######plot change point#######

p1 <-solea_tot[solea_tot$year < 2010,]
mean1 <- mean(p1$AbunRecruits)
p2 <-solea_tot[solea_tot$year>= 2010,]
mean2 <- mean(p2$AbunRecruits)

p <- ggplot()
p <-p+geom_segment(aes(x = 2007, xend = 2010, y = 124, yend = 124), col="black")#create segment
p <-p+geom_segment(aes(x = 2010, xend = 2019, y = 263, yend = 263), col="black")
p <- p+ geom_point(data = solea_tot, aes(x =year, y =AbunRecruits ),color = "seagreen3",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =AbunRecruits ),color = "seagreen3")
p <- p+ ylab("Recruits (thousands)")+xlab("Year")
p <-p+geom_vline(aes(xintercept =2010), col="black", linetype = "dashed")#vertical line
p <- p+ theme_bw()
p
ggsave("RecruitsCP_sole.png", p, width=10, height=10, dpi=300,units="cm" )


p1 <-solea_tot[solea_tot$year < 2013,]
mean1 <- mean(p1$AbunAdults)
p2 <-solea_tot[solea_tot$year>= 2013,]
mean2 <- mean(p2$AbunAdults)

p <- ggplot()
p <-p+geom_segment(aes(x = 2007, xend = 2013, y = 159, yend = 159), col="black")#create segment
p <-p+geom_segment(aes(x = 2013, xend = 2019, y = 390, yend = 390), col="black")
#p <- p+ geom_line(data = data_sol_y, aes(x =year2, y =new_w_sqkm ),color = "black",size=1)
p <- p+ geom_point(data = solea_tot, aes(x =year, y =AbunAdults ),color = "#41b6c4",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =AbunAdults ),color = "#41b6c4")
p <-p+geom_vline(aes(xintercept =2013), col="black", linetype = "dashed")#vertical line
p <- p+ ylab("Adults (thousands)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("AdultsCP_sole.png", p, width=10, height=10, dpi=300,units="cm" )

p1 <-solea_tot[solea_tot$year < 2012,]
mean1 <- mean(p1$Abun_tot)
p2 <-solea_tot[solea_tot$year>= 2012,]
mean2 <- mean(p2$Abun_tot)

p <- ggplot()
p <-p+geom_segment(aes(x = 2007, xend = 2012, y = 291, yend = 291), col="black")#create segment
p <-p+geom_segment(aes(x = 2012, xend = 2019, y = 646, yend = 646), col="black")
#p <- p+ geom_line(data = data_sol_y, aes(x =year2, y =new_w_sqkm ),color = "black",size=1)
p <- p+ geom_point(data = solea_tot, aes(x =year, y =Abun_tot ),color = "#225ea8",size=1)
p <- p+geom_smooth(data = solea_tot, aes(x =year, y =Abun_tot ),color = "#225ea8")
p <-p+geom_vline(aes(xintercept =2012), col="black", linetype = "dashed")#vertical line
p <- p+ ylab("Sole (number/sqkm)")+xlab("Year")
p <- p+ theme_bw()
p
ggsave("TotCP_sole.png", p, width=10, height=10, dpi=300,units="cm" )


###########bimodality#######

#version plot1 
a <- ggplot(solea_tot, aes(x = AbunRecruits))
p <- a + geom_density(col="seagreen3") +
        geom_vline(aes(xintercept = mean(AbunRecruits)), 
             linetype = "dashed", size = 0.6,col="seagreen3")+
  xlab("Recruits (thousands)")+ylab("Density")+
  theme_bw()
p
ggsave("bim_r.png", p, width=10, height=10, dpi=300,units="cm" )

a <- ggplot(solea_tot, aes(x = AbunAdults))
p <- a + geom_density(color = "#41b6c4") +
  geom_vline(aes(xintercept = mean(AbunAdults)), 
             linetype = "dashed", size = 0.6,color = "#41b6c4")+
  xlab("Adults (thousands)")+ylab("Density")+
  theme_bw()
p
ggsave("bim_adult.png", p, width=10, height=10, dpi=300,units="cm" )

a <- ggplot(solea_tot, aes(x = Abun_tot))
p <- a + geom_density(color = "#225ea8") +
  geom_vline(aes(xintercept = mean(Abun_tot)), 
             linetype = "dashed", size = 0.6,color = "#225ea8")+
  xlab("Total Sole (number/sqkm)")+ylab("Density")+
  theme_bw()
p
ggsave("bim_tot.png", p, width=10, height=10, dpi=300,units="cm" )


####### driver state plot########

########state -driver plot 

solea_tot$cp_r <- NA
solea_tot$cp_r[solea_tot$year<2010]<-1
solea_tot$cp_r[solea_tot$year>=2010]<-2

solea_tot$cp_r <- as.factor(solea_tot$cp_r)

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =bottom_t, y =AbunRecruits),col="seagreen3")
p <- p+ geom_point(data = solea_tot, aes(x =bottom_t, y =AbunRecruits, col=cp_r),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Sea Bottom Temperature (°C)")+ylab("Recruits (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_rec_SBT.png", p, width=10, height=10, dpi=300,units="cm" )


p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =runoff_ms, y =AbunRecruits),col="seagreen3")
p <- p+ geom_point(data = solea_tot, aes(x =runoff_ms, y =AbunRecruits, col=cp_r),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Runoff Po River")+ylab("Recruits (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_rec_runoff.png", p, width=10, height=10, dpi=300,units="cm" )

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =PP, y =AbunRecruits),col="seagreen3")
p <- p+ geom_point(data = solea_tot, aes(x =PP, y =AbunRecruits, col=cp_r),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Primary Production")+ylab("Recruits (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_rec_PP.png", p, width=10, height=10, dpi=300,units="cm" )

##adults
solea_tot$cp_A <- NA
solea_tot$cp_A[solea_tot$year<2013]<-1
solea_tot$cp_A[solea_tot$year>=2013]<-2

solea_tot$cp_A <- as.factor(solea_tot$cp_A)

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =bottom_t, y =AbunAdults),col="#41b6c4")
p <- p+ geom_point(data = solea_tot, aes(x =bottom_t, y =AbunAdults, col=cp_A),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Sea Bottom Temperature (°C)")+ylab("Adults (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_adults_SBT.png", p, width=10, height=10, dpi=300,units="cm" )


p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =runoff_ms, y =AbunAdults),col="#41b6c4")
p <- p+ geom_point(data = solea_tot, aes(x =runoff_ms, y =AbunAdults, col=cp_A),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Runoff Po River")+ylab("Adults (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_adults_runoff.png", p, width=10, height=10, dpi=300,units="cm" )

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =effort, y =AbunAdults),col="#41b6c4")
p <- p+ geom_point(data = solea_tot, aes(x =effort, y =AbunAdults, col=cp_A),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Effort")+ylab("Adults (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_adults_effort.png", p, width=10, height=10, dpi=300,units="cm" )

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =PP, y =AbunAdults),col="#41b6c4")
p <- p+ geom_point(data = solea_tot, aes(x =PP, y =AbunAdults, col=cp_A),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Primary Production")+ylab("Adults (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_adults_PP.png", p, width=10, height=10, dpi=300,units="cm" )

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =AbunRecruits, y =AbunAdults),col="#41b6c4")
p <- p+ geom_point(data = solea_tot, aes(x =AbunRecruits, y =AbunAdults, col=cp_A),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Recruits (thousands t)")+ylab("Adults (thousands)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_adults_rec.png", p, width=10, height=10, dpi=300,units="cm" )


#tot
solea_tot$cp_T <- NA
solea_tot$cp_T[solea_tot$year<2012]<-1
solea_tot$cp_T[solea_tot$year>=2012]<-2

solea_tot$cp_T <- as.factor(solea_tot$cp_T)

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =bottom_t, y =Abun_tot),col="#225ea8")
p <- p+ geom_point(data = solea_tot, aes(x =bottom_t, y =Abun_tot, col=cp_T),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Sea Bottom Temperature (°C)")+ylab("Total Sole (number/sqkm)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_tot_SBT.png", p, width=10, height=10, dpi=300,units="cm" )


p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =runoff_ms, y =Abun_tot),col="#225ea8")
p <- p+ geom_point(data = solea_tot, aes(x =runoff_ms, y =Abun_tot, col=cp_T),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Runoff Po River")+ylab("Total Sole (number/sqkm)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_tot_runoff.png", p, width=10, height=10, dpi=300,units="cm" )

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =effort, y =Abun_tot),col="#225ea8")
p <- p+ geom_point(data = solea_tot, aes(x =effort, y =Abun_tot, col=cp_T),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Effort")+ylab("Total Sole (number/sqkm)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_tot_effort.png", p, width=10, height=10, dpi=300,units="cm" )

p <- ggplot()
p <- p+geom_path(data = solea_tot, aes(x =PP, y =Abun_tot),col="#225ea8")
p <- p+ geom_point(data = solea_tot, aes(x =PP, y =Abun_tot, col=cp_T),size=3)+scale_color_manual(values = c("#00AFBB",  "#FC4E07"))
p <- p+xlab("Primary Production")+ylab("Total Sole (number/sqkm)")
p <- p+theme_bw()
p <-p + theme(legend.position = "none")
p
ggsave("plot_state_driver_tot_PP.png", p, width=10, height=10, dpi=300,units="cm" )

##########Modelling the total abundance########
m1 <- lm(Abun_tot ~ bottom_t, data=solea_tot)
summary(m1)
plot(m1)
#plot(m1)
pred <-predict(m1)
sol <- cbind(solea_tot, pred)
p1 <- ggplot(data = sol,
             aes(bottom_t, Abun_tot)) + geom_point() + 
  geom_smooth(method="lm", col="#225ea8", se=T)+theme_bw()+ xlab("Sea Bottom Temperature (°C)")+
  ylab("Total Sole (number/sqkm)")
p1
ggsave("model_SBT_tot.png", p1, width=10, height=10, dpi=300,units="cm" )

#threshold gam

library(mgcv)#gam package
library(INDperform)#tGAM package
library(readxl)
library(readr)
library(tidyverse)
library(patchwork)


y <- solea_tot$Abun_tot

x <- solea_tot$runoff_ms

x2 <- solea_tot$effort

time <- solea_tot$year

mod <- gam(y ~ s(x,k=3)) 

tmod <- thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                   k = 4, a = 0.2, b = 0.8)                             
#test_interaction
loocv_thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                 k = 4, a = 0.2, b = 0.8,time =  time ) #[1] FALSE

summary(tmod)

tmod$mr # 9.789377

tmod$train_na <- "NA"
tmod$train_na <- rep(FALSE, times = 56) 

plot_diagnostics(tmod)$all_plots

solea_tot$CP_E <- ifelse(solea_tot$effort>=33735839, 1, 2)
solea_max <- solea_tot[solea_tot$CP_E == 1,]
solea_min <- solea_tot[solea_tot$CP_E ==2,]
solea_tot$col <- ifelse(solea_tot$CP_E==1, "#00AFBB",  "#FC4E07")

p <- ggplot(data = solea_tot, aes(x = runoff_ms, y = Abun_tot)) + 
  geom_point(color=solea_tot$col) +
  geom_smooth(data=solea_max,aes(x = runoff_ms, y = Abun_tot), method = "lm", se = TRUE, col="#00AFBB")+
  geom_smooth(data=solea_min,aes(x = runoff_ms, y = Abun_tot),method = "lm", se = TRUE, col="#FC4E07")+
  ylab("Total Sole (number/sqkm)")+xlab("Runoff Po River")+
  theme_bw()
p

ggsave("model_runoff_tot_effort.png", p, width=10, height=10, dpi=300,units="cm" )


##########Modelling Recruitment#######

m1 <- lm(AbunRecruits ~ bottom_t, data=solea_tot)
summary(m1)
plot(m1)
pred <-predict(m1)
p1 <- ggplot(data = cbind(solea_tot, pred),
             aes(bottom_t, AbunRecruits)) + geom_point() + geom_line(aes(y=pred))+theme_bw()+ xlab("SBT (C)")+
  ylab("Recruits abundance")
p1
ggsave("model_SBT_R.png", p1, width=10, height=10, dpi=300,units="cm" )
m1 <- lm(AbunRecruits ~ runoff_ms, data=solea_tot)
summary(m1)
#no direct effect
m1b <- lm(AbunRecruits ~runoff_ms+cp_r, data=solea_tot)
summary(m1b)
m1a <- lm(AbunRecruits ~bottom_t+cp_r, data=solea_tot)
summary(m1a)
m1c <- lm(AbunRecruits ~runoff_ms*cp_r, data=solea_tot)
summary(m1a)
m1d <- lm(AbunRecruits ~bottom_t*cp_r, data=solea_tot)
summary(m1d)
m2 <- lm(AbunRecruits ~PP, data=solea_tot)
summary(m2)

dat<- solea_tot[,c(1,2,7)]

dat$lag <- c(NA, dat$AbunAdults[1:12])

m3 <- lm(AbunRecruits ~lag, data=dat)
summary(m3)

m1 <- gam(AbunRecruits ~ s(PP,k=3), data=solea_tot)
summary(m1)
plot(m1)
m1 <- gam(AbunRecruits ~ s(runoff_ms,k=3), data=solea_tot)
summary(m1)
plot(m1)

m1 <- gam(AbunRecruits ~ s(bottom_t,k=3), data=solea_tot)
summary(m1)
plot(m1)


######Modelling adults####
m1 <- lm(AbunAdults ~ bottom_t, data=solea_tot)
summary(m1)
plot(m1)
pred <-predict(m1)
p1 <- ggplot(data = cbind(solea_tot, pred),
             aes(bottom_t, AbunAdults)) + geom_point() + geom_smooth(method="lm", col="#41b6c4", se=T)+theme_bw()+ xlab("Sea Bottom Temperature (°C)")+
  ylab("Adults (thousands)")
p1
ggsave("model_SBT_adults.png", p1, width=10, height=10, dpi=300,units="cm" )

#########tgam########

library(mgcv)#gam package
library(INDperform)#tGAM package
library(readxl)
library(readr)
library(tidyverse)
library(patchwork)

solea_tot$abunA_lag <- c(NA,solea_tot$AbunAdults[1:12] )
solea_tot2 <- solea_tot[!is.na(solea_tot$abunA_lag),]

y <- solea_tot2$abunA_lag

x <- solea_tot2$runoff_ms

x2 <- solea_tot2$AbunRecruits

time <- solea_tot2$year

mod <- gam(y ~ s(x,k=3)) 

tmod <- thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                   k = 4, a = 0.2, b = 0.8)                             
#test_interaction
loocv_thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                 k = 4, a = 0.2, b = 0.8,time =  time ) #[1] FALSE

summary(tmod)

tmod$mr # 9.789377

tmod$train_na <- "NA"
tmod$train_na <- rep(FALSE, times = 56) 

plot_diagnostics(tmod)$all_plots



solea_tot$CP_REC <-ifelse(solea_tot$AbunRecruits <240.3324, 1, 2)
solea_tot$col <- ifelse(solea_tot$CP_REC== 1, "#00AFBB",  "#FC4E07")
solea_min <- solea_tot[solea_tot$CP_REC == 1,]
solea_max<- solea_tot[solea_tot$CP_REC ==2,]

p <- ggplot(data = solea_tot, aes(x = runoff_ms, y = AbunAdults)) + 
  geom_point(color=solea_tot$col) +
  geom_smooth(data=solea_max,aes(x = runoff_ms, y = AbunAdults), method = "lm",se = T, col="#FC4E07")+
  geom_smooth(data=solea_min,aes(x = runoff_ms, y = AbunAdults), method = "lm", se = T, col="#00AFBB")+
  ylab("Adults (thousands)")+xlab("Runoff Po River")+
  theme_bw()
p

ggsave("model_runoff_abundanceadults_runoff_recruits.png", p, width=10, height=10, dpi=300,units="cm" )

y <- solea_tot$AbunAdults

x <- solea_tot$runoff_ms

x2 <- solea_tot$effort

time <- solea_tot$year

mod <- gam(y ~ s(x,k=3)) 

tmod <- thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                   k = 4, a = 0.2, b = 0.8)                             
#test_interaction
loocv_thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                 k = 4, a = 0.2, b = 0.8,time =  time ) #[1] FALSE

summary(tmod)

tmod$mr # 9.789377

tmod$train_na <- "NA"
tmod$train_na <- rep(FALSE, times = 56) 

plot_diagnostics(tmod)$all_plots



solea_tot$CP_EFF_A <-ifelse(solea_tot$effort <33735839, 2, 1)
solea_tot$col <- ifelse(solea_tot$CP_EFF_A== 1, "#00AFBB",  "#FC4E07")
solea_min <- solea_tot2[solea_tot$CP_EFF_A == 1,]
solea_max<- solea_tot2[solea_tot$CP_EFF_A ==2,]

p <- ggplot(data = solea_tot, aes(x = runoff_ms, y = AbunAdults)) + 
  geom_point(color=solea_tot$col) +
  geom_smooth(data=solea_max,aes(x = runoff_ms, y = AbunAdults), method = "lm", se = T, col="#FC4E07")+
  geom_smooth(data=solea_min,aes(x = runoff_ms, y = AbunAdults), method = "lm", se = T, col="#00AFBB")+
  ylab("Adults (thousands)")+xlab("Runoff Po River")+
  theme_bw()
p

ggsave("model_runoff_abundanceadults_runoff_recruits_lag.png", p, width=10, height=10, dpi=300,units="cm" )


y <- solea_tot$AbunRecruits

x <- solea_tot$PP

x2 <- solea_tot$bottom_t

time <- solea_tot$year

mod <- gam(y ~ s(x,k=3)) 

tmod <- thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                   k = 4, a = 0.2, b = 0.8)                             
#test_interaction
loocv_thresh_gam(model = mod, ind_vec = y, press_vec = x, t_var = x2, name_t_var = "x2",
                 k = 4, a = 0.2, b = 0.8,time =  time ) #[1] FALSE

summary(tmod)

tmod$mr # 9.789377

tmod$train_na <- "NA"
tmod$train_na <- rep(FALSE, times = 56) 

plot_diagnostics(tmod)$all_plots



solea_tot$CP_EFF_A <-ifelse(solea_tot$effort <33735839, 2, 1)
solea_tot$col <- ifelse(solea_tot$CP_EFF_A== 1, "#00AFBB",  "#FC4E07")
solea_min <- solea_tot2[solea_tot$CP_EFF_A == 1,]
solea_max<- solea_tot2[solea_tot$CP_EFF_A ==2,]

p <- ggplot(data = solea_tot, aes(x = runoff_ms, y = AbunAdults)) + 
  geom_point(color=solea_tot$col) +
  geom_smooth(data=solea_max,aes(x = runoff_ms, y = AbunAdults), method = "lm", se = FALSE, col="#FC4E07")+
  geom_smooth(data=solea_min,aes(x = runoff_ms, y = AbunAdults), method = "lm", se = FALSE, col="#00AFBB")+
  ylab("Adults (thousands)")+xlab("Runoff Po River")+
  theme_bw()
p

ggsave("model_runoff_abundanceadults_runoff_effort_lag.png", p, width=10, height=10, dpi=300,units="cm" )

