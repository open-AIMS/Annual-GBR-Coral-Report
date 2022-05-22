library(ggsn) #map features
source('CoralTrends_functions.R')
CoralTrends_checkPackages()
source('CoralTrends_config.R')

load('../data/processed/manta.sum.RData')
load(file='../data/modelled/dat.gbr.RData')
#load(file='../data/modelled/last_year.RData')
load(file='../data/modelled/mod.gbr.RData')

load(file='../data/spatial/whagbr.RData')
load(file='../data/spatial/whagbr.n.RData')
load(file='../data/spatial/whagbr.c.RData')
load(file='../data/spatial/whagbr.s.RData')
load(file='../data/spatial/qld.RData')

## Before looking at modelled data, lets make a plot of raw data for the finalyear

ewbrks <- seq(144,152,by=2)
nsbrks <- seq(-24,-10,by=2)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste0(x, "°W"), ifelse(x > 0, paste0(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), "°S"), ifelse(x > 0, paste0(x, "°N"),x))))
trafficLightColors = colorRampPalette(c('#ED1C24','#F47721','#F0C918','#B0D235','#00734D'))

g = manta.sum %>% filter(REPORT_YEAR==final_year) %>% droplevels %>%
  mutate(cCover = cut(Cover, breaks=c(0,0.10,0.30,0.50,0.75,1), labels=1:5)) %>%
  ggplot() +
  geom_polygon(data=fortify(qld), aes(y=lat, x=long,group=group), fill='grey', color='grey40') +
  geom_polygon(data=fortify(whagbr.n),aes(y=lat, x=long, group=group), color='black', fill=NA) +
  geom_polygon(data=fortify(whagbr.c),aes(y=lat, x=long, group=group), color='black', fill=NA) +
  geom_polygon(data=fortify(whagbr.s),aes(y=lat, x=long, group=group), color='black', fill=NA) +
  #geom_point(aes(y=Latitude,x=Longitude,fill=cCover, size=Tows),alpha=1,shape=21,color='black', show.legend = TRUE)  +
  geom_point(aes(y=Latitude,x=Longitude,fill=cCover),size=3, alpha=1,shape=21,color='black', show.legend = TRUE)  +
  scale_fill_manual('Hard coral cover', breaks=1:5, values=trafficLightColors(5), labels=c('0-10%','10-30%','30-50%','50-75%','75-100%'),limits=1:5) + 
  coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
  annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
  scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
  scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))
  annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
  annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
  annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
  theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                                        #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                          legend.position = c(0.975,0.65),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA)) +
  ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
  north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12) +
  guides(fill=guide_legend(override.aes = list(size=5)))
g
ggsave(file='../output/figures/CoralCoverBubbleMap_OnlyLastYearRawReefs.png', g, width=7, height=7, dpi=300)
ggsave(file='../output/figures/CoralCoverBubbleMap_OnlyLastYearRawReefs.pdf', g, width=7, height=7, dpi=300)
manta.sum %>% filter(REPORT_YEAR==final_year) %>% droplevels %>%
    mutate(cCover = cut(Cover, breaks=c(0,0.1,0.3,0.5,0.75,1), labels=1:5)) %>%
    write_csv(file='../data/processed/FinalYear_coral_cover.csv')


## Now we can look at the change between the final_year and the previous observation (per reef)
last_year = manta.sum %>% group_by(REEF_NAME) %>%
      summarize(MaxYr=max(REPORT_YEAR, na.rm=TRUE),
                MaxYrCover=Cover[REPORT_YEAR==MaxYr],
                MinYr=max(REPORT_YEAR[REPORT_YEAR<MaxYr], na.rm=TRUE),
                MinYrCover=Cover[REPORT_YEAR==MinYr],
                Longitude=mean(Longitude, na.rm=TRUE),
                Latitude=mean(Latitude, na.rm=TRUE)) %>%
      filter(MinYr>=(final_year-2)) %>%
      mutate(DiffYr=MaxYr-MinYr,
             Diff = (MaxYrCover - MinYrCover)/DiffYr,
             D=Diff<0)

ewbrks <- seq(144,152,by=2)
nsbrks <- seq(-24,-10,by=2)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste0(x, "°W"), ifelse(x > 0, paste0(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), "°S"), ifelse(x > 0, paste0(x, "°N"),x))))
manta.sum.reefs = manta.sum %>% filter(REPORT_YEAR==final_year) %>% dplyr:::select(REEF_NAME) %>% distinct

ly=last_year %>% right_join(manta.sum.reefs) %>% filter(!is.na(MaxYr))
save(ly, file='../data/modelled/ly.RData')
ly %>% write.csv(file='../data/processed/FinalYear_coral_change.csv')

gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    ##    geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    ##    geom_point(data=ly, aes(y=Latitude,x=Longitude,size=abs(Diff)*100), shape=21, alpha=0, color='black') +
    ##geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    ##geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    ##geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    #geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=D), size=abs(ly$Diff)*ly$Tows,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
                                        #geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=D), size=abs(ly$Diff)*100/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    ##    geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=D, size=abs(Diff)*100),alpha=0.5,shape=21,color='black', show.legend = TRUE) +
    geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=D, size=abs(Diff)*100),alpha=0.5,shape=21,color='black', show.legend = TRUE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('green','red'),limits=c(FALSE,TRUE))+
    ##scale_size('Magnitude of change', breaks=c(4,8,12,16), labels=c(4,8,12,16),values=c(4,8,12,16))+
                                        #scale_size_identity('Magnitude of change') +
    scale_size('Magnitude', breaks=c(0.5,1,2,5,10,15),range=c(0.5,15)/2) +
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.975,0.65),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
gp
ggsave(file='../output/figures/ChangeInCoralCoverBubbleMap_OnlyLastYearReefs.png', gp, width=7, height=7, dpi=300)
ggsave(file='../output/figures/ChangeInCoralCoverBubbleMap_OnlyLastYearReefs.pdf', gp, width=7, height=7, dpi=300)


## For 2020, Mike also wanted the annual change year by year ==================================================================
## This is based on a figure made in the production of Mikes paper
## Note, this differs from above in that it does not restrict the time span to only two years back....
load('../data/processed/manta.sum.RData')
management <- whagbr.n + whagbr.c + whagbr.s
management.df <- fortify(management)
tops = management.df %>% filter(lat > max(lat)-0.1) %>% summarize(min=min(long),max=max(long), lat=max(lat))
manta.change=manta.sum  %>%
  mutate(Location=Region) %>% 
  group_by(Location,Zone,REEF_NAME,Longitude,Latitude) %>%
  arrange(REPORT_YEAR) %>%
  mutate(Cover.lag = lag(Cover),
         Change = 100*(Cover-Cover.lag),#/Cover.lag,
         Time.lag=lag(REPORT_YEAR),
         Time.change=REPORT_YEAR-Time.lag,
         Change.annual=Change/Time.change) %>%
    filter(REPORT_YEAR>min(REPORT_YEAR)) %>%
  ungroup %>%
  droplevels
gA=ggplot(data=NULL, aes(y=Latitude,x=Longitude)) #+
    #geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2)
long.offset=8
lat.offset=-17
Yrs = sort(unique(manta.change$REPORT_YEAR))
NCOL=7
NROW=ceiling(length(Yrs)/NCOL)
m=1:(NCOL*NROW)
m=m[m>NCOL*length(Yrs)/NCOL]
mat=matrix(c(sort(unique(manta.change$REPORT_YEAR))), ncol=NCOL, byrow=TRUE)
tmat = t(mat)
tmat[m] <- NA
mat = t(tmat)
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    if (is.na(mat[i,j])) next
    gA=gA +
      geom_polygon(data=management.df %>% mutate(lat=lat+(i-1)*lat.offset, long=long+(j-1)*long.offset), aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
      geom_point(data=manta.change %>% filter(REPORT_YEAR==mat[i,j]) %>% mutate(Latitude=Latitude+(i-1)*lat.offset, Longitude=Longitude+(j-1)*long.offset),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
      geom_text(data=tops %>% mutate(lat=lat+0.2 + (i-1)*lat.offset, min=min+(j-1)*long.offset, max=max+(j-1)*long.offset, long=mean(c(min,max))), aes(y=lat, x=long),label=mat[i,j], vjust=0) 
  }
}

gA=gA + coord_equal() +
    scale_fill_manual('',breaks=c(FALSE,TRUE), values=c('red','green'), labels=c('Decrease','Increase')) +
    scale_radius('Annual\nchange (%)', breaks=c(5,10,15,20,25)) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.line=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    guides(fill=guide_legend(order=1, label.theme=element_text(size=10), override.aes=list(size=3)),
           size=guide_legend(order=0, label.theme=element_text(size=10)))
gA
ggsave(filename=paste0('../output/figures/AnnualChange_',3,'Zones_new.pdf'), gA, width=10,height=10.2)
save(manta.change, file='../data/processed/manta.change.RData')


## For 2020, Mike also wants a bubble plot of COTS and Bleaching
## COTS
## The following data set was provided by Mike via email on 29/06/2020
## NOTE, THIS IS TEMPORARY FOR 2020, THIS SHOULD BE FULLY CODED FROM NOW ON..
## NOTE IN PARTICULAR, THAT WADE REEF IS MANUALLY CODED HERE TO ACCOUNT FOR A NA THAT SHOULD NOT HAVE BEEN THERE.
##cots <- read.csv('../data/primary/COTS per reef-AIMS-8HTC3Z2.csv', strip.white=TRUE)
## For 2021, Mike has supplied both COTS and Bleaching in one file.
## This file uses different mappings of COTS categories
cotsAndBleaching <- read_csv('../data/primary/COTS and bleaching for annual report.csv')
## Mike supplied the following conversions
##COTS = 0: No COTS
##COTS>0 & COTS<0.1: No Outbreak
##COTS>=0.1 & COTS<0.22: Outbreak watch
##COTS>=0.22 & COTS<=1: Incipient Outbreak
##COTS>1: Active Outbreak
##In addition, if STATUS==RE: Recovering
cots <- cotsAndBleaching %>%
  mutate(OUTBREAK.CAT = case_when(
           MEAN_COTS==0 ~ 'No COTS',
           MEAN_COTS>0 & MEAN_COTS<0.1  ~ 'No Outbreak',
           MEAN_COTS>=0.1 & MEAN_COTS<0.22 ~ 'Outbreak Watch',
           MEAN_COTS>=0.22 & MEAN_COTS<=1 ~ 'Incipient Outbreak',
           MEAN_COTS>1 ~ 'Active Outbreak',
           STATUS=='RE' ~ 'Recovering'
         )
         ) %>%
  mutate(OUTBREAK.CAT = ifelse(OUTBREAK.CAT=='Outbreak watch', 'Outbreak Watch', OUTBREAK.CAT),
         ## OUTBREAK.CAT = ifelse(REEF_NAME=='WADE REEF', 'Outbreak Watch', OUTBREAK.CAT),
         OUTBREAK.CAT = factor(OUTBREAK.CAT, levels=c('No COTS', 'No Outbreak','Outbreak Watch', 'Recovering', 'Incipient Outbreak', 'Active Outbreak')))
## Actually, it turns out that 'Outbreak Watch' should be 'Potential Outbreak'
cots = cots %>%
  mutate(OUTBREAK.CAT = as.character(OUTBREAK.CAT),
         OUTBREAK.CAT = ifelse(OUTBREAK.CAT=='Outbreak Watch', 'Potential Outbreak', OUTBREAK.CAT),
         OUTBREAK.CAT = factor(OUTBREAK.CAT, levels=c('No COTS', 'No Outbreak','Potential Outbreak', 'Recovering', 'Incipient Outbreak', 'Active Outbreak')))

gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=cots, aes(y=REEF_LAT,x=REEF_LONG,fill=OUTBREAK.CAT),alpha=1,shape=21,size=3,color='black', show.legend = TRUE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('COTS', breaks=levels(cots$OUTBREAK.CAT), labels=levels(cots$OUTBREAK.CAT), values=rev(trafficLightColors(6)),limits=levels(cots$OUTBREAK.CAT))+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.975,0.65),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=1)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
gp
ggsave(file='../output/figures/COTSBubbleMap_OnlyLastYearReefs.png', gp, width=7, height=7, dpi=300)
ggsave(file='../output/figures/COTSBubbleMap_OnlyLastYearReefs.pdf', gp, width=7, height=7, dpi=300)

## Bleaching
## bleaching <- readxl::read_excel('../data/primary/aethstetic bleaching.xlsx')
bleaching <- cotsAndBleaching %>%
  mutate(`Bleaching Aes` = MaxOfBLEACHING_PERHC)
bleaching <- bleaching %>%
  mutate(CAT=ifelse(is.na(`Bleaching Aes`), 0,
             ifelse(`Bleaching Aes`==0, 0,
             ifelse(`Bleaching Aes` %in% c('0+', '1-', '1', '1+'), '<10%',
             ifelse(`Bleaching Aes` %in% c('2-', '2', '2+'), '10-30%',
             ifelse(`Bleaching Aes` %in% c('3-', '3', '3+'), '30-50%','>50%')
             )))),
         ## CAT = factor(CAT, levels=c(NA,'0','<10%','10-30%','30-50%','>50%'))) %>%
         CAT = factor(CAT, levels=c('0','<10%','10-30%','30-50%','>50%'))) %>%
  left_join(cots %>% dplyr::select(REEF_LAT, REEF_LONG, REEF_NAME) %>% distinct)

gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
  geom_polygon(aes(group=group), fill='grey', color='grey40') +
  geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
  geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
  geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
  geom_point(data=bleaching, aes(y=REEF_LAT,x=REEF_LONG,fill=CAT),alpha=1,shape=21,size=3,color='black', show.legend = TRUE) +
  coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
  annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
  scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
  scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
  scale_fill_manual('Bleaching',
                    breaks=c(levels(bleaching$CAT)), #breaks=c(NA,levels(bleaching$CAT)),
                    labels=c(levels(bleaching$CAT)), #labels=c('NA',levels(bleaching$CAT)),
                    values=c(rev(trafficLightColors(5))), #values=c('grey',rev(trafficLightColors(5))),
                    limits=c(levels(bleaching$CAT))) + #limits=c(NA,levels(bleaching$CAT)))+
  annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
  annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
  annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
  theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                                        #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                          legend.position = c(0.975,0.65),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
  guides(fill=guide_legend(override.aes = list(size=5,alpha=1)),
         size=guide_legend(override.aes = list(alpha=0.5)),
         color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                       dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
  north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
gp


## Mike would like a four panel plot with current coral cover, change in cover, cots and bleaching ========================================================================
## To do this, we must:
## 1. make the base figure (qld, sites etc)
gbr.sp<- readShapeSpatial("~/Work/Resources/GIS/Features/Great_Barrier_Reef_Features.shp",
                          proj4string = CRS('+proj=longlat +ellps=WGS84'),repair=TRUE,force_ring=T,verbose=TRUE)
Zones=c('Northern','Central','Southern')
## Clipped coast
bb = bbox(management)
##gbr.sp = raster:::crop(gbr.sp, raster::extent(bb))
gbr.fort <- broom::tidy(gbr.sp, region='FEAT_NAME')
grd.y = data.frame(x=c(140,140,140,140,140,140,140),
                 xend=c(146,146,148,150,152,154,154),
                 y=c(-11,-13,-15,-17,-19,-21,-23),
                 yend=c(-11,-13,-15,-17,-19,-21,-23))
grd.x = data.frame(x=c(142,144,146,148,150,152,154),
                 xend=c(142,144,146,148,150,152,154),
                 y=c(-9,-9,-11,-15,-17,-19,-21),
                 yend=c(-30,-30,-30,-30,-30,-30,-30))

g.base = ggplot(manta.sum %>% filter(REPORT_YEAR==final_year), aes(y=Latitude,x=Longitude)) +
  geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long, group=group),fill='grey', color='grey70') +
  geom_polygon(data=fortify(qld), aes(y=lat, x=long, group=group),fill='grey', color='grey40') +
  geom_segment(data=grd.y, aes(x=x,y=y,xend=xend,yend=yend), size=0.1, color='grey80')+
  geom_segment(data=grd.x, aes(x=x,y=y,xend=xend,yend=yend), size=0.1, color='grey80') +
  annotate(geom='text', x=-Inf,y=seq(-23,-11,by=2), label=degreeLabelsNS(seq(-23,-11,by=2)), hjust=-0.1,parse=TRUE, size=2) +
  annotate(geom='text', y=-Inf,x=seq(142,154,by=2), label=degreeLabelsEW(seq(142,154,by=2)),vjust=-0.2,parse=TRUE,size=2) +
  geom_point(shape=21,fill='yellow',color='black') +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.5) +
  geom_point(data=towns12 %>% filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')), aes(y=lat, x=long)) +
  geom_text(data=towns12 %>% filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')), aes(y=lat, x=long, label=town), hjust=1.1, size=3) +
  annotate(geom='text',y=-13,x=145.2,label=Zones[1], angle=-76, hjust=0.5,vjust=-0.5) +
  annotate(geom='text',y=-18,x=147.8,label=Zones[2], angle=-31, hjust=0.5,vjust=-0.5) +
  annotate(geom='text',y=-22,x=153.2,label=Zones[3], angle=-73, hjust=0.5,vjust=-0.5) +
  scale_x_continuous(breaks=seq(142,154,by=2), position = 'bottom', expand=c(0,0)) +
  scale_y_continuous(breaks=seq(-25,-10,by=2)) +
  coord_equal(ylim=bbox(management)[2,],xlim=c(bbox(management)[1,1],bbox(management)[1,2]+20.5)) +
  theme_minimal() +
  theme(legend.position=c(0.05,0.10), legend.justification=c(0,0),
        legend.background = element_rect(color='black', fill='white'),
        axis.title=element_blank(),
        panel.border = element_rect(color='black',fill=NA),
        axis.text=element_blank(),
                                        #axis.ticks.length = unit(c(-0.5,-0.9),'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat=management.df%>% filter(lat< -15 & lat> -16, long>145.5) #%>% filter(lat==first(lat) | lat==last(lat))
(coefs=coef(lm(long~lat,dat)))
strt = c(145.9,-14.5)
ff=function(strt,coefs,length) {
    x=strt[1]
    y=strt[2]
    dat=data.frame(x=x,y=y)
    for (i in 2:length) {
        y=y-0.7
        x=coefs[1] + y*coefs[2]
        dat=rbind(dat, data.frame(x=x[[1]], y=y[[1]]))
    }
    dat
}

## 2. Coral cover for the current year
lgnd.dat = ff(c(coefs[1] + (strt[2]*coefs[2]),strt[2]),coefs,length=5)
lgnd.dat$lab = rev(c('>0% - 10%','>10% - 30%','>30% - 50%','>50% - 75%','>75% - 100%'))
lgnd.dat$Cat=factor(5:1)

g.cover = manta.sum %>%
  filter(REPORT_YEAR==final_year) %>% droplevels %>%
  mutate(cCover = cut(Cover, breaks=c(0,0.10,0.30,0.50,0.75,1), labels=1:5)) %>%
  ggplot(aes(y=Latitude, x=Longitude)) +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
  geom_point(aes(y=Latitude,x=Longitude,fill=cCover),size=3, alpha=1,shape=21,color='black', show.legend = TRUE)  +
  scale_fill_manual('Hard coral cover', breaks=1:5, values=trafficLightColors(5), labels=rev(c('>0% - 10%','>10% - 30%','>30% - 50%','>50% - 75%','>75% - 100%')),limits=1:5) +
  geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='Coral cover'), vjust=0) +
  scale_x_continuous(expand=c(0,0)) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.4, fill=Cat), size=3, shape=21) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 
g.cover
g.cover=ggplotGrob(g.cover)
g.cover.grob.panel=g.cover[[1]][[6]]
shiftacross = 5
g = g.base + annotation_custom(grob=g.cover.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 

## Change in coral cover (between last year and the previous within 3 years)
lgnd.dat = data.frame(x=c(145.2,145.2),y=c(-12,-12.7)) #ff(strt,coefs,length=2)
lgnd.dat$lab = c("Increase","Decrease")
lgnd.dat$Cat=factor(c(FALSE,TRUE))
## lgnd.dat.size = ff(lgnd.dat[2,1:2],coefs,length=5)
lgnd.dat.size = ff(c(coefs[1] + (strt[2]*coefs[2]), strt[2]),coefs,length=5)
lgnd.dat.size$lab = format(seq(2.5,12.5,by=2.5),digits=2)
lgnd.dat.size$size = seq(2.5,12.5,by=2.5)

g.change = ly %>%
  ggplot(aes(y=Latitude, x=Longitude)) +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
  geom_point(aes(y=Latitude,x=Longitude,fill=D, size=abs(Diff*100)),alpha=1,shape=21,color='black', show.legend = TRUE)  +
  scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('green','red'),limits=c(FALSE,TRUE)) +
  geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='Coral change'), vjust=0) +
  scale_x_continuous(expand=c(0,0)) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=TRUE) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_point(data=lgnd.dat.size, aes(y=y,x=x+0.4, size=size), shape=21) +
  geom_text(data=lgnd.dat.size, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0, parse=TRUE) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())    +
  scale_radius()
g.change
g.change=ggplotGrob(g.change)
g.change.grob.panel=g.change[[1]][[6]]
shiftacross = shiftacross + 5
g = g + annotation_custom(grob=g.change.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 

## COTS
##strt = c(145.9,-14.5)
lgnd.dat = ff(c(coefs[1] + (-14.5+0.7)*coefs[2], -14.5+0.7),coefs,length=6)
lgnd.dat$lab = levels(cots$OUTBREAK.CAT)
lgnd.dat$Cat=factor(lgnd.dat$lab)
g.cots = cots %>%
  dplyr::rename(Latitude=REEF_LAT, Longitude=REEF_LONG)%>%
  ggplot(aes(y=Latitude, x=Longitude)) +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
  geom_point(aes(y=Latitude,x=Longitude,fill=OUTBREAK.CAT),size=3, alpha=1,shape=21,color='black', show.legend = TRUE) +
  scale_fill_manual('COTS', breaks=levels(cots$OUTBREAK.CAT), labels=levels(cots$OUTBREAK.CAT), values=rev(trafficLightColors(6)),limits=levels(cots$OUTBREAK.CAT)) +
  geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='COTS'), vjust=0) +
  scale_x_continuous(expand=c(0,0)) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=FALSE) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
g.cots
g.cots=ggplotGrob(g.cots)
g.cots.grob.panel=g.cots[[1]][[6]]
shiftacross = shiftacross + 5
g = g + annotation_custom(grob=g.cots.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 

## Bleaching
levels(bleaching$CAT) <-  c('0%', '>0% - 10%', '>10% - 30%', '>30% - 50%', '>50%') 
lgnd.dat = ff(c(coefs[1] + (strt[2]*coefs[2]), strt[2]),coefs,length=5)
lgnd.dat$lab = levels(bleaching$CAT)
lgnd.dat$Cat=factor(lgnd.dat$lab)
g.bleaching = bleaching %>%
  dplyr::rename(Latitude=REEF_LAT, Longitude=REEF_LONG)%>%
  ggplot(aes(y=Latitude, x=Longitude)) +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
  geom_point(aes(y=Latitude,x=Longitude,fill=CAT),size=3, alpha=1,shape=21,color='black', show.legend = TRUE) +
  scale_fill_manual('Bleaching',
                    breaks=c(NA,levels(bleaching$CAT)),
                    labels=c('NA',levels(bleaching$CAT)),
                    values=c('grey',rev(trafficLightColors(5))),
                    limits=c(NA,levels(bleaching$CAT))) +
  geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='Bleaching'), vjust=0) +
  scale_x_continuous(expand=c(0,0)) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=FALSE) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
g.bleaching
g.bleaching=ggplotGrob(g.bleaching)
g.bleaching.grob.panel=g.bleaching[[1]][[6]]
shiftacross = shiftacross + 5
g = g + annotation_custom(grob=g.bleaching.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 

## Add the Australia inset
coord.df <- map_data("world2", "Australia", exact=FALSE,
                  xlim=c(110,155),ylim=c(-45,-5),
                  boundary=FALSE,
                  interior=TRUE, as.polygon=TRUE)
coord <- maps:::map.poly("world", "Australia", exact=FALSE,
                         xlim=c(110,160),ylim=c(-45,-5),
                         boundary=FALSE,
                         interior=TRUE, fill=TRUE, as.polygon=TRUE)
IDs <- sapply(strsplit(coord$names, ":"), function(x) x[1])
coord.sp <- map2SpatialPolygons(coord,IDs=IDs)
bb.4=bb
#bb.4[1,1] = bb[1,1]-(bb[1,2] - bb[1,1])*0.1
bb.4[2,1] = bb[2,1]-(bb[2,2] - bb[2,1])*0.1
gClip <- function(shp, bb, type='Intersection'){
    require(raster)
    require(rgeos)
    if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
    else b_poly <- as(extent(bb), "SpatialPolygons")
    if(type=='Intersection') return(gIntersection(shp, b_poly, byid = T))
    if(type=='Difference') return(gDifference(b_poly,shp, byid = T))
}

coord.sp = gClip(coord.sp, bb.4)

aus=ggplot(coord.df, aes(y=lat,x=long,group=group)) + geom_polygon(fill='white',color='black') + coord_map() +
    geom_polygon(data=fortify(coord.sp), aes(y=lat, x=long, group=group),fill='grey',color='black', size=0.2) +
    geom_polygon(data=fortify(management), aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
                                        #annotate(geom='rect', xmin=bb[1,1],xmax=bb[1,2], ymin=bb[2,1],ymax=bb[2,2], fill=NA, color='black') +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          panel.background = element_rect(color='black'))
aus.grob=ggplotGrob(aus)
g = g + annotation_custom(grob=aus.grob,xmax=Inf,ymax=Inf,xmin=169,ymin=-15)


## Add the scalebar and north arrow
gc.1=c(143.5,-24)
gc.2=gcDestination(lon=gc.1[1], lat=gc.1[2],bearing=90, dist=500, model='WGS84', Vincenty=FALSE)
g = g + ggsn:::scalebar(x.min=gc.1[1],y.min=gc.1[2],x.max=gc.2[1],y.max=gc.1[2]+0.5,
                        dist=250, model='WGS84',st.size=2, height=0.5,st.dist=0.5, dist_unit = 'km',
                        transform=TRUE) +
    ggsn:::north(x.min=mean(c(gc.1[1],gc.2[1]))-0.5, x.max=mean(c(gc.1[1],gc.2[1]))+0.5, y.min=gc.1[2]+0.5,y.max=gc.1[2]+2, scale=1) 
g1=g
g1


library(magick)
library(png)
a=magick::image_read(path='../parameters/ECM_1280725_v1_AIMS Logo - stacked.jpg') 

g1 <- g1 + annotation_custom(rasterGrob(a, x=unit(0.92, 'npc'), y=unit(0.52,'npc'), width=unit(0.1, 'npc')))

ggsave(filename=paste0('../output/figures/FourPlots.pdf'), g1, width=12,height=5.8)
ggsave(filename=paste0('../output/figures/FourPlots.png'), g1, width=12,height=6.8, dpi=300)

## Mike also wants a version of this figure without the sampling locations - this will reduce the width of the figure ======================
## He also wants to have the sampling dates appended to the 3 region names
tops = management.df %>% filter(lat > max(lat)-0.1) %>% summarize(min=min(long),max=max(long), lat=max(lat))
lgnd.dat = ff(c(coefs[1] + (strt[2]*coefs[2]), strt[2]),coefs,length=5)
lgnd.dat$lab = rev(c('>0% - 10%','>10% - 30%','>30% - 50%','>50% - 75%','>75% - 100%'))
lgnd.dat$Cat=factor(5:1)


#Zones1 = paste0(Zones, c('\n(Nov-Dec 2019)', '\n(Feb/June 2020)', '\n(Sept/Oct 2019 & Jan 2020)'))
Zones1 = c('Survey dates\n(Oct-Dec 2020)', 'Survey dates\n(Jan/Feb/April 2021)', 'Survey dates\n(Aug 2020)')

g.base = manta.sum %>% filter(REPORT_YEAR==final_year) %>%
  mutate(cCover = cut(Cover, breaks=c(0,0.10,0.30,0.50,0.75,1), labels=1:5)) %>%
  ggplot(, aes(y=Latitude,x=Longitude)) +
  geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long, group=group),fill='grey', color='grey70') +
  geom_polygon(data=fortify(qld), aes(y=lat, x=long, group=group),fill='grey', color='grey40') +
  geom_segment(data=grd.y, aes(x=x,y=y,xend=xend,yend=yend), size=0.1, color='grey80')+
  geom_segment(data=grd.x, aes(x=x,y=y,xend=xend,yend=yend), size=0.1, color='grey80') +
  annotate(geom='text', x=-Inf,y=seq(-23,-11,by=2), label=degreeLabelsNS(seq(-23,-11,by=2)), hjust=-0.1,parse=TRUE, size=2) +
  annotate(geom='text', y=-Inf,x=seq(142,154,by=2), label=degreeLabelsEW(seq(142,154,by=2)),vjust=-0.2,parse=TRUE,size=2) +
  #geom_point(shape=21,fill='yellow',color='black') +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.5) +
  geom_point(aes(y=Latitude,x=Longitude,fill=cCover),size=3, alpha=1,shape=21,color='black', show.legend = TRUE)  +
  scale_fill_manual('Hard coral cover', breaks=1:5, values=trafficLightColors(5), labels=rev(c('>0% - 10%','>10% - 30%','>30% - 50%','>50% - 75%','>75% - 100%')),limits=1:5) +
  annotate(geom='rect', xmin=tops$min, xmax=tops$max, ymin=tops$lat+0.01, ymax=tops$lat+0.7, fill='white', alpha=0.8) + 
  geom_text(data=tops %>% mutate(min=min, max=max, long=min,max), aes(y=lat+0.1, x=long+0.1,label='a) Coral cover'), vjust=0, hjust=0) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.4, fill=Cat), size=3, shape=21) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0) +
  geom_point(data=towns12 %>% filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')), aes(y=lat, x=long)) +
  geom_text(data=towns12 %>% filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')), aes(y=lat, x=long, label=town), hjust=1.1, size=3) +
                                        #annotate(geom='text',y=-12.8,x=145.2,label=Zones1[1], angle=-76, hjust=0.5,vjust=-0.5, size=2.5) +
                                        #annotate(geom='text',y=-19,x=149.4,label=Zones1[2], angle=-31, hjust=0.5,vjust=-0.5, size=2.5) +
                                        #annotate(geom='text',y=-22.8,x=153.4,label=Zones1[3], angle=-73, hjust=0.5,vjust=-0.5, size=2.5) +
  annotate(geom='text',y=-13,x=145.2,label=Zones[1], angle=-76, hjust=0.5,vjust=-0.5) +
  annotate(geom='text',y=-18.5,x=148.7,label=Zones[2], angle=-31, hjust=0.5,vjust=-0.5) +
  annotate(geom='text',y=-22,x=153.3,label=Zones[3], angle=-73, hjust=0.5,vjust=-0.5) +
  scale_x_continuous(breaks=seq(142,154,by=2), position = 'bottom', expand=c(0,0)) +
  scale_y_continuous(breaks=seq(-25,-10,by=2)) +
  coord_equal(ylim=bbox(management)[2,],xlim=c(bbox(management)[1,1],bbox(management)[1,2]+16)) +
  theme_minimal() +
  theme(legend.position='none', #legend.position=c(0.05,0.10), legend.justification=c(0,0),
        legend.background = element_rect(color='black', fill='white'),
        axis.title=element_blank(),
        panel.border = element_rect(color='black',fill=NA),
        axis.text=element_blank(),
                                        #axis.ticks.length = unit(c(-0.5,-0.9),'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g.base

# Coral change
lgnd.dat = data.frame(x=c(145.2,145.2),y=c(-12,-12.7)) #ff(strt,coefs,length=2)
lgnd.dat$lab = c("Increase","Decrease")
lgnd.dat$Cat=factor(c(FALSE,TRUE))
## lgnd.dat.size = ff(lgnd.dat[2,1:2],coefs,length=5)
lgnd.dat.size = ff(c(coefs[1] + (strt[2]*coefs[2]), strt[2]),coefs,length=5)
lgnd.dat.size$lab = format(seq(2.5,12.5,by=2.5),digits=2)
lgnd.dat.size$size = seq(2.5,12.5,by=2.5)

g.change = ly %>%
  ggplot(aes(y=Latitude, x=Longitude)) +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
  geom_point(aes(y=Latitude,x=Longitude,fill=D, size=abs(Diff*100)),alpha=1,shape=21,color='black', show.legend = TRUE)  +
  scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('green','red'),limits=c(FALSE,TRUE)) +
  geom_text(data=tops %>% mutate(min=min, max=max, long=min,max), aes(y=lat+0.1, x=long+0.1,label='b) Coral change'), vjust=0, hjust=0) +
  scale_x_continuous(expand=c(0,0)) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=TRUE) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_point(data=lgnd.dat.size, aes(y=y,x=x+0.4, size=size), shape=21) +
  geom_text(data=lgnd.dat.size, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0, parse=TRUE) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())    +
  scale_radius()
g.change=ggplotGrob(g.change)
g.change.grob.panel=g.change[[1]][[6]]
shiftacross = 5
g = g.base + annotation_custom(grob=g.change.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 

## COTS
lgnd.dat = ff(c(coefs[1] + (strt[2]+0.7)*coefs[2], strt[2]+0.7),coefs,length=6)
lgnd.dat$lab = levels(cots$OUTBREAK.CAT)
lgnd.dat$Cat=factor(lgnd.dat$lab)
g.cots = cots %>%
  dplyr::rename(Latitude=REEF_LAT, Longitude=REEF_LONG)%>%
  ggplot(aes(y=Latitude, x=Longitude)) +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
  geom_point(aes(y=Latitude,x=Longitude,fill=OUTBREAK.CAT),size=3, alpha=1,shape=21,color='black', show.legend = TRUE) +
  scale_fill_manual('COTS', breaks=levels(cots$OUTBREAK.CAT), labels=levels(cots$OUTBREAK.CAT), values=rev(trafficLightColors(6)),limits=levels(cots$OUTBREAK.CAT)) +
  geom_text(data=tops %>% mutate(min=min, max=max, long=min,max), aes(y=lat+0.1, x=long+0.1,label='c) COTS'), vjust=0, hjust=0) +
  scale_x_continuous(expand=c(0,0)) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=FALSE) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
g.cots
g.cots=ggplotGrob(g.cots)
g.cots.grob.panel=g.cots[[1]][[6]]
shiftacross = shiftacross + 5
g = g + annotation_custom(grob=g.cots.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 

## Bleaching
levels(bleaching$CAT) <-  c('0%', '>0% - 10%', '>10% - 30%', '>30% - 50%', '>50%') 
lgnd.dat = ff(c(coefs[1] + strt[2]*coefs[2], strt[2]),coefs,length=5)
lgnd.dat$lab = levels(bleaching$CAT)
lgnd.dat$Cat=factor(lgnd.dat$lab)
g.bleaching = bleaching %>%
  dplyr::rename(Latitude=REEF_LAT, Longitude=REEF_LONG)%>%
  ggplot(aes(y=Latitude, x=Longitude)) +
  geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
  geom_point(aes(y=Latitude,x=Longitude,fill=CAT),size=3, alpha=1,shape=21,color='black', show.legend = TRUE) +
  annotate(geom='text',y=-12.8,x=145.2,label=Zones1[1], angle=-76, hjust=0.5,vjust=-0.5, size=3) +
  annotate(geom='text',y=-19,x=149.4,label=Zones1[2], angle=-31, hjust=0.5,vjust=-0.5, size=3) +
  annotate(geom='text',y=-22.8,x=153.4,label=Zones1[3], angle=-73, hjust=0.5,vjust=-0.5, size=3) +
  scale_fill_manual('Bleaching',
                    breaks=c(levels(bleaching$CAT)),
                    labels=c(levels(bleaching$CAT)),
                    values=c(rev(trafficLightColors(5))),
                    limits=c(levels(bleaching$CAT))) +
  geom_text(data=tops %>% mutate(min=min, max=max, long=min,max), aes(y=lat+0.1, x=long+0.1,label='d) Bleaching'), vjust=0, hjust=0) +
  scale_x_continuous(expand=c(0,0)) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=FALSE) +
  geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
g.bleaching
g.bleaching=ggplotGrob(g.bleaching)
g.bleaching.grob.panel=g.bleaching[[1]][[6]]
shiftacross = shiftacross + 5
g = g + annotation_custom(grob=g.bleaching.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 

## Add the Australia inset
coord.df <- map_data("world2", "Australia", exact=FALSE,
                  xlim=c(110,155),ylim=c(-45,-5),
                  boundary=FALSE,
                  interior=TRUE, as.polygon=TRUE)
coord <- maps:::map.poly("world", "Australia", exact=FALSE,
                         xlim=c(110,160),ylim=c(-45,-5),
                         boundary=FALSE,
                         interior=TRUE, fill=TRUE, as.polygon=TRUE)
IDs <- sapply(strsplit(coord$names, ":"), function(x) x[1])
coord.sp <- map2SpatialPolygons(coord,IDs=IDs)
bb.4=bb
#bb.4[1,1] = bb[1,1]-(bb[1,2] - bb[1,1])*0.1
bb.4[2,1] = bb[2,1]-(bb[2,2] - bb[2,1])*0.1
gClip <- function(shp, bb, type='Intersection'){
    require(raster)
    require(rgeos)
    if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
    else b_poly <- as(extent(bb), "SpatialPolygons")
    if(type=='Intersection') return(gIntersection(shp, b_poly, byid = T))
    if(type=='Difference') return(gDifference(b_poly,shp, byid = T))
}

coord.sp = gClip(coord.sp, bb.4)

aus=ggplot(coord.df, aes(y=lat,x=long,group=group)) + geom_polygon(fill='white',color='black') + coord_map() +
    geom_polygon(data=fortify(coord.sp), aes(y=lat, x=long, group=group),fill='grey',color='black', size=0.2) +
    geom_polygon(data=fortify(management), aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
                                        #annotate(geom='rect', xmin=bb[1,1],xmax=bb[1,2], ymin=bb[2,1],ymax=bb[2,2], fill=NA, color='black') +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          panel.background = element_rect(color='black'))
aus.grob=ggplotGrob(aus)
g = g + annotation_custom(grob=aus.grob,xmax=Inf,ymax=Inf,xmin=164,ymin=-15)


## Add the scalebar and north arrow
gc.1=c(143.5,-24)
gc.2=gcDestination(lon=gc.1[1], lat=gc.1[2],bearing=90, dist=500, model='WGS84', Vincenty=FALSE)
g = g + ggsn:::scalebar(x.min=gc.1[1],y.min=gc.1[2],x.max=gc.2[1],y.max=gc.1[2]+0.5,
                        dist=250, model='WGS84',st.size=2, height=0.5,st.dist=0.5, dist_unit = 'km',
                        transform=TRUE) +
    ggsn:::north(x.min=mean(c(gc.1[1],gc.2[1]))-0.5, x.max=mean(c(gc.1[1],gc.2[1]))+0.5, y.min=gc.1[2]+0.5,y.max=gc.1[2]+2, scale=1) 
g1=g
g1

library(magick)
library(png)
a=magick::image_read(path='../parameters/ECM_1280725_v1_AIMS Logo - stacked.jpg') 

g1 <- g1 + annotation_custom(rasterGrob(a, x=unit(0.92, 'npc'), y=unit(0.55,'npc'), width=unit(0.1, 'npc')))

ggsave(filename=paste0('../output/figures/FourPlots1.pdf'), g1, width=10,height=5.6)
ggsave(filename=paste0('../output/figures/FourPlots1.png'), g1, width=10,height=6.8, dpi=300)



## STOP HERE (maybe) ===========================================

