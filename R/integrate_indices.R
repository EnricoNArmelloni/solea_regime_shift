library(readxl)


# recruits
rec_index=read_excel("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/data/area_index_below20cm.xlsx")
rec_2009=read_excel("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/data/area_index_below20cm_y2009.xlsx")
rec_index[rec_index$Survey== 'SOLEMON2009_b', 11:14]=rec_2009[,11:14]

# adults
adu_2009=read_excel("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/data/area_index_above20cm_y2009.xlsx")
adu_index=read_excel("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/data/area_index_above20cm.xlsx")
adu_index[adu_index$Survey== 'SOLEMON2009_b', 11:14]=adu_2009[,11:14]


# all
all_index=read_excel("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/data/area_index_all.xlsx")


# combine
all=all_index[,c('Survey','Species','AbunIndex','AbunCV','BiomIndex','BiomCV')]
rec=rec_index[,c('Survey','Species','AbunIndex','AbunCV')]
names(rec)[3:4]=c('AbunRecruits','AbunRecruitsCV')
adu=adu_index[,c('Survey','Species','AbunIndex','AbunCV')]
names(adu)[3:4]=c('AbunAdults','AbunAdultsCV')

combined_index=left_join(all, rec)%>%
  left_join(adu)

combined_index$year=as.numeric(str_remove(str_remove(str_remove(str_remove(combined_index$Survey,'SOLEMON'),'_b'),'OTT'),'NOVEMBRE'))

combined_index=combined_index[combined_index$year<=2019,]

ggplot(data=combined_index)+
  geom_line(aes(x=year, y=AbunIndex))

p1=ggplot(data=combined_index,aes(x=year, y=AbunAdults))+
  geom_line()+
  geom_point()
p2=ggplot(data=combined_index,aes(x=year, y=AbunRecruits))+
  geom_line()+
  geom_point()
p3=ggplot(data=combined_index,aes(x=year, y=AbunIndex))+
  geom_line()+
  geom_point()
p4=ggplot(data=combined_index,aes(x=year, y=BiomIndex))+
  geom_line()+
  geom_point()
ggpubr::ggarrange(p1,p2,p3,p4,ncol=1)

write.csv(combined_index, 'C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/data/Survey_indices.csv',row.names = F)

ggsave('C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/solea_regime_shift/plots/Indices.png', width = 10, height = 20, units='cm')
