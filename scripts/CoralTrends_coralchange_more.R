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


## rstanarm have made a change to ranef.  as of June 2019, it now calls lmer/glmer to use as a template
## from which the values are then replaced with the stanfit
## unfortunately, glmer takes forever with these fits, so I am going to use the old
## methodology
ranef.stanreg <- function(object, ...) {
  all_names <- if (rstanarm:::used.optimizing(object))
    rownames(object$stan_summary) else object$stanfit@sim$fnames_oi
  sel <- rstanarm:::b_names(all_names)
  ans <- object$stan_summary[sel, rstanarm:::select_median(object$algorithm)]
  # avoid returning the extra levels that were included
  ans <- ans[!grepl("_NEW_", names(ans), fixed = TRUE)]
  fl <- rstanarm:::.flist.stanreg(object)
  levs <- lapply(fl, levels)
  asgn <- attr(fl, "assign")
  cnms <- rstanarm:::.cnms.stanreg(object)
  fl <- fl
  asgn <- asgn
  levs <- levs
  cnms <- cnms
  nc <- vapply(cnms, length, 1L)
  nb <- nc * vapply(levs, length, 1L)
  nbseq <- rep.int(seq_along(nb), nb)
  ml <- split(ans, nbseq)
  for (i in seq_along(ml)) {
    ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE, 
                      dimnames = list(NULL, cnms[[i]]))
  }
  ans <- lapply(seq_along(fl), function(i) {
    data.frame(do.call(cbind, ml[i]), row.names = levs[[i]], 
               check.names = FALSE)
  })
  names(ans) <- names(fl)
  structure(ans, class = "ranef.mer")
}

## From here down are modelled estimates and these are not used in the >=2019 figures
CoralTrends_changeInYears = function(mod,yr1,yr2) {
    load('data/processed/manta.sum.RData')    
    a=fixef(mod)[1]
    b=fixef(mod)[-1]
    ## rstanarm have made a change to ranef - see comment above
    a.reef=ranef.stanreg(mod)[[1]][,'(Intercept)']
    b.reef=as.matrix(ranef(mod)[[1]][,-1])
    
    bb=sweep(b.reef,2,b,FUN='+')
    bbb = sweep(bb,1,(a+a.reef), FUN='+')
    state = data.frame(binomial()$linkinv((bbb)))
    state$REEF_NAME = rownames(state)
    state = state %>% left_join(manta.sum %>% group_by(REEF_NAME) %>% summarize_at(vars(Latitude,Longitude,Tows), funs(mean)))

    ##Determine the most recent year and the next most recent year and keep only
    ##those within two years of the current year
    manta.yrs=manta.sum %>% group_by(REEF_NAME) %>%
      summarize(MaxYr=max(REPORT_YEAR), MinYr=max(REPORT_YEAR[REPORT_YEAR<MaxYr])) %>%
      filter(MinYr>=(final_year-2)) %>%
      mutate(DiffYr=MaxYr-MinYr)

    last_year = state %>% left_join(manta.yrs) %>% filter(!is.na(MaxYr)) %>% group_by(REEF_NAME) %>%
      mutate_(.dots=setNames(paste0('Year',.$MaxYr,'-Year',.$MinYr),'Diff')) %>% 
      dplyr::select(REEF_NAME,Diff,DiffYr,Latitude,Longitude,Tows) %>%
      mutate(Diff=Diff/DiffYr) %>% 
      #mutate_(.dots=setNames(paste0('Year',.$MaxYr,'-Year',.$MaxYr),'Diff')) %>%
      mutate(D=Diff<0) #%>% dplyr:::select(-starts_with('Year'))
    
    ## last_year = state %>% dplyr:::select_('REEF_NAME',paste0('Year',yr2),paste0('Year',yr1)) %>%
    ##     mutate_(.dots=setNames(paste0('Year',yr2,'-Year',yr1),'Diff')) %>%
    ##     mutate(D=Diff<0) %>% dplyr:::select(-starts_with('Year'))
    ## last_year = last_year %>% left_join(manta.sum %>% group_by(REEF_NAME) %>% summarize_at(vars(Latitude,Longitude,Tows), funs(mean)))
    last_year
}

CoralTrends_predict_reef_year = function(mod) {
    load('data/processed/manta.sum.RData')    
    a=fixef(mod)[1]
    b=fixef(mod)[-1]
    a.reef=ranef.stanreg(mod)[[1]][,'(Intercept)']
    b.reef=as.matrix(ranef(mod)[[1]][,-1])
    
    bb=sweep(b.reef,2,b,FUN='+')
    bbb = sweep(bb,1,(a+a.reef), FUN='+')
    state = data.frame(binomial()$linkinv((bbb)))
    state$REEF_NAME = rownames(state)
    state = state %>% left_join(manta.sum %>% group_by(REEF_NAME) %>% summarize_at(vars(Latitude,Longitude,Tows), funs(mean)))
    state
}

## In the past, we had region wide settings for the years to use for the most
## recent change.  For example, for the northern region, mod.northern,yr1 was
## 2015 and mod.northern,yr2 was 2017.  These settings were defined in
## parameters/CoralTrends.conf.
## This year (2019), Mike would like to have it compare the most recent survey
## to the next most recent survey (provided that is not earlier than 2 years
## from present.

load(file='../data/modelled/mod.northern.RData')
last_year.northern = CoralTrends_changeInYears(mod.northern,yr1=coralchange_northern.y1,yr2=coralchange_northern.y2)

load(file='../data/modelled/mod.central.RData')
last_year.central = CoralTrends_changeInYears(mod.central,yr1=coralchange_central.y1,yr2=coralchange_central.y2)

load(file='../data/modelled/mod.southern.RData')
last_year.southern = CoralTrends_changeInYears(mod.southern,yr1=coralchange_southern.y1,yr2=coralchange_southern.y2)

last_year = rbind(last_year.northern, last_year.central, last_year.southern)

ewbrks <- seq(144,152,by=2)
nsbrks <- seq(-24,-10,by=2)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste0(x, "°W"), ifelse(x > 0, paste0(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), "°S"), ifelse(x > 0, paste0(x, "°N"),x))))

ly=last_year %>% right_join(manta.sum.reefs)
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
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.png', gp, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.pdf', gp, width=7, height=7, dpi=300)
data4bubbleplot = ly
save(data4bubbleplot, file='../data/modelled/data4bubbleplot.RData')


gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,size=abs(last_year$Diff)*last_year$Tows), shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    geom_point(data=last_year %>% arrange(desc(D)), aes(y=Latitude,x=Longitude,fill=D), size=2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('blue','red'),limits=c(FALSE,TRUE))+
    #scale_size_area('Influence', max_size=4)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.75),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2], location='topright',scale=0.1,symbol=12)
gp
ggsave(file='../output/figures/InfluenceMap_manta.png', gp, width=5, height=5, dpi=300)
ggsave(file='../output/figures/InfluenceMap_manta.pdf', gp, width=5, height=5, dpi=300)


gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,size=abs(last_year$Diff)*last_year$Tows), shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=abs(last_year$Diff)*last_year$Tows,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('blue','red'),limits=c(FALSE,TRUE))+
    #scale_size_area('Influence', max_size=4)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.75),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2], location='topright',scale=0.1,symbol=12)
gp
ggsave(file='../output/figures/InfluenceMap_mantaBubble.png', gp, width=5, height=5, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble.pdf', gp, width=5, height=5, dpi=300)


## Now we will repeat this, yet only plot the reefs that were actually sampled in the last year
manta.sum.reefs = manta.sum %>% filter(REPORT_YEAR==final_year) %>% dplyr:::select(REEF_NAME) %>% distinct
gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=last_year %>% right_join(manta.sum.reefs), aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,size=abs(last_year$Diff)*last_year$Tows), shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    geom_point(data=last_year %>% right_join(manta.sum.reefs) %>% arrange(desc(D)), aes(y=Latitude,x=Longitude,fill=D), size=2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('blue','red'),limits=c(FALSE,TRUE))+
    #scale_size_area('Influence', max_size=4)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.75),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
  north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2], location='topright',scale=0.1,symbol=12)
gp
ggsave(file='../output/figures/InfluenceMap_manta_OnlyLastYearReefs.png', gp, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_manta_OnlyLastYearReefs.pdf', gp, width=7, height=7, dpi=300)


## The following is the plot used in 2019
ly=last_year %>% right_join(manta.sum.reefs)
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
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.png', gp, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.pdf', gp, width=7, height=7, dpi=300)

## For the 2018 version, they would like the Northern section to indicate Not sampled
gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=D, size=abs(Diff)*100),alpha=0.5,shape=21,color='black', show.legend = TRUE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('green','red'),limits=c(FALSE,TRUE))+
    ##scale_size('Magnitude of change', breaks=c(4,8,12,16), labels=c(4,8,12,16),values=c(4,8,12,16))+
                                        #scale_size_identity('Magnitude of change') +
    scale_size('Magnitude', breaks=c(0.5,1,2,5,10,15),range=c(0.5,15)/2) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.975,0.65),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
## normal years - without the need to blank out Northern region.
choice0=gp + geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
  geom_polygon(aes(group=group), fill='grey', color='grey40') +
  ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
  north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
choice0
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_2017.png', choice0, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_2017.pdf', choice0, width=7, height=7, dpi=300)

## Choice 1
choice1=gp + geom_polygon(data=fortify(whagbr.n),aes(group=group), color='grey90', fill=NA) +
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR\n(not surveyed in 2017/2018)', hjust=0) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
  ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
choice1
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice1.png', choice1, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice1.pdf', choice1, width=7, height=7, dpi=300)

## Choice 2 - this is the choice they have gone for
choice2=gp + geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill='grey90') +
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR\n(not surveyed in 2017/2018)', hjust=0) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
  ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
choice2
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice2.png', choice2, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice2.pdf', choice2, width=7, height=7, dpi=300)

## Choice 3
choice3=gp + geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR\n(not surveyed in 2017/2018)', hjust=0) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
  ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
choice3
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice3.png', choice3, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice3.pdf', choice3, width=7, height=7, dpi=300)
## Choice 4
choice4=gp + geom_polygon(data=fortify(whagbr.n),aes(group=group), color='grey70', fill=NA) +
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR\n(not surveyed in 2017/2018)', hjust=0,color='grey70') +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
  ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2]-0.5, y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2]+0.2, location='topright',scale=0.1,symbol=12)
choice4
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice4.png', choice4, width=7, height=7, dpi=300)
ggsave(file='../output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs_choice4.pdf', choice4, width=7, height=7, dpi=300)


## Now lets do this on the raw data only
## ---- CoralChangeBubble
load('../data/processed/manta.sum.RData')
## Generate a list of reefs we are using to help accumulate other
## sources of associated data
all.reefs = manta.sum %>%
    dplyr:::select(P_CODE.mod,REEF_NAME,REEF_ID,Latitude,Longitude) %>%
    group_by(REEF_NAME,REEF_ID) %>%
    summarize_at(vars(Latitude,Longitude), funs(mean)) %>%
    as.data.frame

## Genuine stan cannot handle proportional data for binomial families
## (particularly when weights are applied). A work-around is to
## multiple the proportion by the weights and convert this into an integer  
dat.all = manta.sum %>%
    dplyr:::select(Cover, REEF_NAME, Tows,P_CODE.mod,Location,REPORT_YEAR,Latitude,Longitude) %>%
    mutate(Year=factor(REPORT_YEAR), N=length(unique(REEF_NAME))) %>% ungroup() %>%
    mutate(Cvr1 = as.integer(as.vector(Cover) * Tows), Cvr0 = Tows - Cvr1)

## Calculate some relatively simple comparisons
dat.all.gbr = dat.all %>% droplevels

## Define Pre and Post
## dat.all.gbr = dat.all.gbr %>%
##     mutate(Comp='None',
##            Comp = ifelse(Location=='Northern' & REPORT_YEAR==2013, 'Pre', Comp),
##            Comp = ifelse(Location=='Northern' & REPORT_YEAR==2017, 'Post', Comp),
##            Comp = ifelse(Location=='Central' & REPORT_YEAR==2016, 'Pre', Comp),
##            Comp = ifelse(Location=='Central' & REPORT_YEAR==2018, 'Post', Comp),
##            Comp = ifelse(Location=='Southern' & REPORT_YEAR==2016, 'Pre', Comp),
##            Comp = ifelse(Location=='Southern' & REPORT_YEAR==2018, 'Post', Comp)
##            ) %>% filter(Comp!='None') %>%
##     mutate(Location=factor(Location, levels=c('Northern','Central','Southern')))
dat.all.gbr = dat.all.gbr %>%
    mutate(Comp='None',
           Comp = ifelse(Location=='Northern' & REPORT_YEAR==2017, 'Pre', Comp),
           Comp = ifelse(Location=='Northern' & REPORT_YEAR==2019, 'Post', Comp),
           Comp = ifelse(Location=='Central' & REPORT_YEAR==2017, 'Pre', Comp),
           Comp = ifelse(Location=='Central' & REPORT_YEAR==2019, 'Post', Comp),
           Comp = ifelse(Location=='Southern' & REPORT_YEAR==2017, 'Pre', Comp),
           Comp = ifelse(Location=='Southern' & REPORT_YEAR==2019, 'Post', Comp)
           ) %>% filter(Comp!='None') %>%
    mutate(Location=factor(Location, levels=c('Northern','Central','Southern')))

dd=dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Pre','Post'))) %>%
    dplyr::select(Latitude,Longitude,Tows,Cover,Comp)

gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=dd,
               aes(y=Latitude,x=Longitude,size=abs(Tows), fill=Cover),
               alpha=0.5,shape=21,color='black', show.legend = TRUE) +
    facet_grid(~Comp) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.65),
                            legend.justification=c(1,0.5),
                            legend.background = element_rect(fill='white', color=NA))
    
gp


gp1=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=dd %>% filter(Comp=='Pre'),
               aes(y=Latitude,x=Longitude,size=Tows, fill=Cover*100),
               alpha=0.5,shape=21,color='black', show.legend = TRUE) +
    scale_size_area('Tows', breaks=c(25,50,75,100), limits=c(0,100)) +
    scale_fill_gradient('% Coral Cover') +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.65),
                            legend.justification=c(1,0.5),
                            legend.box='horizontal',
                            legend.background = element_rect(fill='white', color=NA),
                            panel.background=element_rect(fill=NA, color='black')) +
    ggtitle('Pre')
gp1=gp1 + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                   dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
  north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
        location='topright',scale=0.1,symbol=12)


gp2=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=dd %>% filter(Comp=='Post'),
               aes(y=Latitude,x=Longitude,size=Tows, fill=Cover*100),
               alpha=0.5,shape=21,color='black', show.legend = TRUE) +
    scale_size_area('Tows', breaks=c(25,50,75,100), limits=c(0,100)) +
    scale_fill_gradient('% Coral Cover') +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.65),
                            legend.justification=c(1,0.5),
                            legend.box='horizontal',
                            legend.background = element_rect(fill='white', color=NA),
                            panel.background=element_rect(fill=NA, color='black')) +
    ggtitle('Post')
gp2=gp2 + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                   dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
          location='topright',scale=0.1,symbol=12)
 
dd.1=dat.all.gbr %>% group_by(Location, Latitude,Longitude,REEF_NAME) %>%
     mutate(Tows=median(Tows,na.rm=TRUE)) %>%
     dplyr::select(Location,Latitude,Longitude,REEF_NAME,Tows,Comp,Cover) %>%
     spread(key=Comp, value=Cover) %>%
     ungroup %>%
     mutate(Diff=Pre-Post, DiffP=100*Diff/Pre) %>%
     dplyr::select(Latitude,Longitude,Tows,DiffP) %>% filter(!is.na(DiffP))

gp3=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=dd.1,
               aes(y=Latitude,x=Longitude,size=Tows, fill=DiffP*-1),
               alpha=0.5,shape=21,color='black', show.legend = TRUE) +
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_gradient2('% Change',low='red',mid='white',high='green')+
    scale_size_area('Tows', breaks=c(25,50,75,100), limits=c(0,100)) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.65),
                            legend.justification=c(1,0.5),
                            legend.box='horizontal',
                            legend.background = element_rect(fill='white', color=NA),
                            panel.background=element_rect(fill=NA, color='black')) +
    ggtitle('Pre vs Post')
gp3=gp3 + ggsn::scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                   dist=200, model = 'WGS84',st.size=3,location='bottomleft', transform=TRUE, dist_unit = 'km') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
          location='topright',scale=0.1,symbol=12)


library(gridExtra)
grid.arrange(gp1, gp2, gp3, nrow=1)
ggsave(file='../output/figures/CoralChange_bubble.png', grid.arrange(gp1, gp2, gp3, nrow=1), width=15, height=6, dpi=300)
ggsave(file='../output/figures/CoralChange_bubble.pdf', grid.arrange(gp1, gp2, gp3, nrow=1), width=15, height=6, dpi=300)
## ----
