source('CoralTrends_functions.R') 
CoralTrends_checkPackages()
source('CoralTrends_config.R')

INCLUDE_GBR <- TRUE

## ---- loadData
load('../data/processed/manta.sum.RData')
load(file='../data/modelled/dat.gbr.RData')
load(file='../data/modelled/dat.northern.RData')
load(file='../data/modelled/dat.central.RData')
load(file='../data/modelled/dat.southern.RData')

load(file='../data/spatial/spatial_3Zone.RData')

#load(file='data/modelled/cots.sum.all_3Zone.RData')
#load(file='data/modelled/bleaching.sum.all_3Zone.RData')
#load(file='data/modelled/cyclones.sum.all_3Zone.RData')
## ----end


###################################################################################
## The plots will have tick marks along the temporal (x) axis.                   ##
## For the forseable future, these should represent every five years.            ##
## The starting year will always be 1985, but the final year will keep evolving. ##
## Lets take the variable final_year (defined in parameters/CoralTrends.conf     ##
## and round it to the nearest 5 and use that.                                   ##
###################################################################################
mceiling <- function(x,base){ 
        base*ceiling(x/base) 
} 
final_year_seq <- mceiling(final_year,5)

#Start with all panel plot
newdata <- rbind(dat.gbr,dat.northern,dat.central,dat.southern)
newdata$Location <- factor(newdata$Location, levels=unique(newdata$Location),
                           labels=c('Great Barrier Reef','Northern GBR', 'Central GBR', 'Southern GBR'))

if (!INCLUDE_GBR) newdata <- newdata %>% filter(Location!='Great Barrier\n\nReef') %>% droplevels
#newdata = newdata %>% left_join(dat.all %>% select(Location, N) %>% distinct)
nd <- newdata %>%
    group_by(Location) %>%
    summarize(Year=mean(range(as.numeric(as.character(Year)))),
              N=paste0('(N=',unique(N),")"))
hues <- RColorBrewer::brewer.pal(4, "Blues")

g1<-ggplot(newdata, aes(y = mean*100, x = as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow = 1, scales='fixed', labeller = labeller(Location = setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2])+
    #geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x = as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,60)) +
    scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
    theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color=NA),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
                                        #panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                                        #panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
          strip.text=element_text(margin=margin(t=1, b=1,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.7,vjust=-1),
          plot.margin=unit(c(0,0,2,0),'pt'))
write.csv(newdata, file='output/Modelled_coral_cover.csv',quote=FALSE, row.names=FALSE)
          
## Now we generate the banner
a=oz:::ozRegion(sections=c(3,11:13))
a=oz:::ozRegion()
cc=rbind(xy2df(a$lines[[3]]),
    xy2df(a$lines[[13]]),
    xy2df(a$lines[[12]])[nrow(xy2df(a$lines[[12]])):1,],
    xy2df(a$lines[[11]]))
aa.ps<-SpatialPolygons(list(Polygons(list(Polygon(cc)),ID="QLD")))

gt = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
    #geom_blank(aes(x=190,y=-20))+coord_map() +
    #geom_blank(aes(x=220,y=-20))+coord_map() +
    geom_blank(aes(x=270,y=-20))+coord_map() +
    geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
    coord_equal() +
    theme_classic() +theme(panel.background=element_rect(fill=NA),
                       axis.text.y=element_blank(),
                       axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                       axis.ticks=element_blank(),
                       axis.line=element_blank(),
                       plot.background=element_blank(),
                       panel.spacing=unit(0,'pt'),
                       plot.margin=unit(c(0,0,0,0),'pt'))  

gt1=gt+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4]) #+
    #annotate(geom='text', x=Inf, y=Inf, label='Far\nNorthern', vjust=1,hjust=1)

gt2=gt+geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt3=gt+geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt4=gt+geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])


if (!INCLUDE_GBR) {gt1=gt2; gt2=gt3; gt3=gt4;}

## finally we put them together
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=6, name="pic_predator"))
facets <- grep("strip-t-2-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=9, b=b, r=10, name="pic_predator"))
facets <- grep("strip-t-3-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=13, b=b, r=14, name="pic_predator"))
##facets <- grep("strip-t-4-1", gT$layout$name)
##gg <- with(gg$layout[facets,],
##           gtable_add_grob(gg, ggplotGrob(gt4),t=t, l=16, b=b, r=16, name="pic_predator"))
if (INCLUDE_GBR) {
    facets <- grep("strip-t-4-1", gT$layout$name)
                                        #gg <- with(gg$layout[facets,],
                                        #           gtable_add_grob(gg, ggplotGrob(gt5),t=t, l=20, b=b, r=20, name="pic_predator"))
    gg <- with(gg$layout[facets,],
               gtable_add_grob(gg, ggplotGrob(gt4),t=t, l=16, b=b, r=18, name="pic_predator"))
}

grid.draw(gg)
save(gg, file='data/spatial/gg_3Zone.RData')
ggsave(file='output/figures/3ZonesStan.pdf', grid.draw(gg), width=15, height=3, units='in',dpi=300) 
ggsave(file='output/figures/3ZonesStan.png', grid.draw(gg), width=15, height=3, units='in',dpi=300) 


## For talks, Mike needs a three panel (North, Central, Southern) with means and
## errorbars (not ribbon)
newdata1 = newdata %>% filter(Location!='Great Barrier Reef') %>% droplevels() %>%
    mutate(Location=factor(Location, levels=c('Northern GBR', 'Central GBR', 'Southern GBR'),
                           labels=c('Northern', 'Central', 'Southern')))
g1<- newdata1 %>%
    ggplot(aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata1$Location),"\n"), levels(newdata1$Location))))+
    geom_blank()+
    #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2])+
    geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    #geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,60)) +
    scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
    theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color=NA),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
                                        #panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                                        #panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
          strip.text=element_text(margin=margin(t=1, b=1,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.7,vjust=-1),
          plot.margin=unit(c(0,7,2,0),'pt'),
          panel.spacing.x=unit(10,'pt'))
g1
gt.new = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
    #geom_blank(aes(x=190,y=-20))+coord_map() +
    #geom_blank(aes(x=220,y=-20))+coord_map() +
    geom_blank(aes(x=200,y=-20))+coord_map() +
    geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
    coord_equal() +
    theme_classic() +theme(panel.background=element_rect(fill=NA),
                       axis.text.y=element_blank(),
                       axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                       axis.ticks=element_blank(),
                       axis.line=element_blank(),
                       plot.background=element_blank(),
                       panel.spacing=unit(0,'pt'),
                       plot.margin=unit(c(0,0,0,0),'pt'))  

gt1.new=gt.new+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4]) #+
    #annotate(geom='text', x=Inf, y=Inf, label='Far\nNorthern', vjust=1,hjust=1)

gt2.new=gt.new+geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt3.new=gt.new+geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt4.new=gt.new+geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

## finally we put them together
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt2.new),t=t, l=5, b=b, r=6, name="pic_northern"))
facets <- grep("strip-t-2-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt3.new),t=t, l=9, b=b, r=10, name="pic_central"))
facets <- grep("strip-t-3-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt4.new),t=t, l=13, b=b, r=14, name="pic_southern"))

grid.draw(gg)
ggsave(file='output/figures/3ZonesStanErrorBars.pdf', grid.draw(gg), width=11.25, height=3, units='in',dpi=300) 
ggsave(file='output/figures/3ZonesStanErrorBars.png', grid.draw(gg), width=11.25, height=3, units='in',dpi=300) 


## now a three panel but with ribbons ======================================================================
g1<- newdata1 %>%
    ggplot(aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata1$Location),"\n"), levels(newdata1$Location))))+
    geom_blank()+
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2])+
    #geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    #geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,60)) +
    scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
    theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color=NA),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
                                        #panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                                        #panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
          strip.text=element_text(margin=margin(t=1, b=1,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.7,vjust=-1),
          plot.margin=unit(c(0,7,2,0),'pt'),
          panel.spacing.x=unit(10,'pt'))
g1
gt.new = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
    #geom_blank(aes(x=190,y=-20))+coord_map() +
    #geom_blank(aes(x=220,y=-20))+coord_map() +
    geom_blank(aes(x=200,y=-20))+coord_map() +
    geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
    coord_equal() +
    theme_classic() +theme(panel.background=element_rect(fill=NA),
                       axis.text.y=element_blank(),
                       axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                       axis.ticks=element_blank(),
                       axis.line=element_blank(),
                       plot.background=element_blank(),
                       panel.spacing=unit(0,'pt'),
                       plot.margin=unit(c(0,0,0,0),'pt'))  

gt1.new=gt.new+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4]) #+
    #annotate(geom='text', x=Inf, y=Inf, label='Far\nNorthern', vjust=1,hjust=1)

gt2.new=gt.new+geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt3.new=gt.new+geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt4.new=gt.new+geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

## finally we put them together
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt2.new),t=t, l=5, b=b, r=6, name="pic_northern"))
facets <- grep("strip-t-2-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt3.new),t=t, l=9, b=b, r=10, name="pic_central"))
facets <- grep("strip-t-3-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt4.new),t=t, l=13, b=b, r=14, name="pic_southern"))

grid.draw(gg)
ggsave(file='output/figures/3ZonesStanRibbons.pdf', grid.draw(gg), width=11.25, height=3, units='in',dpi=300) 
ggsave(file='output/figures/3ZonesStanRibbons.png', grid.draw(gg), width=11.25, height=3, units='in',dpi=300) 



## There is also a need for separate panels
## GBR
g.gbr<-newdata %>% filter(Location=='Great Barrier Reef') %>% ggplot(aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2])+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd %>% filter(Location=='Great Barrier Reef'), aes(y=50,x=Year, label=N), vjust=1.2) +
        scale_y_continuous(expression(Coral~cover~('%'))) +
            scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
                theme_classic()+
                    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                          panel.background=element_rect(color='black'),
                          plot.margin = margin(t=2,r=7,b=0,l=0),
                          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
                          axis.text.x=element_text(size=rel(1.2)),
                          axis.text.y=element_text(size=rel(1.2)),
                          panel.grid.minor=element_line(size=0.1,color='gray70'),
                          panel.grid.major=element_line(size=0.1,color='gray70'),
                          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
                          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
                          strip.text=element_text(margin=margin(t=1, b=1, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
ggsave(file='output/figures/mantaCAT_gbrNoBanner.png', g.gbr + facet_null(), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_gbrNoBanner.pdf', g.gbr + facet_null(), width=5, height=3.5, units='in',dpi=300)

gt1=gt+geom_polygon(data=fortify(whagbr) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(whagbr), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])
g <- ggplot_gtable(ggplot_build(g.gbr))
facets <- grep("strip-t-1-1", g$layout$name)
gg.gbr <- with(g$layout[facets,],
           gtable_add_grob(g, ggplotGrob(gt1),t=t, l=5, b=b, r=6, name="pic_predator"))

ggsave(file='output/figures/mantaCAT_gbr.png', grid.draw(gg.gbr), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_gbr.pdf', grid.draw(gg.gbr), width=5, height=3.5, units='in',dpi=300)


## Northern GBR
g.northern<-newdata %>% filter(Location=='Northern GBR') %>% ggplot(aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2])+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd %>% filter(Location=='Northern GBR'), aes(y=50,x=Year, label=N), vjust=1.2) +
    scale_y_continuous(expression(Coral~cover~('%'))) +
    scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
    theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
          plot.margin = margin(t=2,r=7,b=0,l=0),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.text.x=element_text(size=rel(1.2)),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color='gray70'),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
          strip.text=element_text(margin=margin(t=1, b=1, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
ggsave(file='output/figures/mantaCAT_northernNoBanner.png', g.northern + facet_null(), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_northernNoBanner.pdf', g.northern + facet_null(), width=5, height=3.5, units='in',dpi=300)

gt2=gt+geom_polygon(data=fortify(whagbr.n), aes(y=lat, x=long),fill=hues[4],color=NA)+
    geom_polygon(data=fortify(whagbr.n), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])
g <- ggplot_gtable(ggplot_build(g.northern))
facets <- grep("strip-t-1-1", g$layout$name)
gg.northern <- with(g$layout[facets,],
                    gtable_add_grob(g, ggplotGrob(gt2),t=t, l=5, b=b, r=6, name="pic_predator"))

ggsave(file='output/figures/mantaCAT_northern.png', grid.draw(gg.northern), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_northern.pdf', grid.draw(gg.northern), width=5, height=3.5, units='in',dpi=300)


## Central GBR
g.central<-newdata %>% filter(Location=='Central GBR') %>% ggplot(aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2])+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd %>% filter(Location=='Central GBR'), aes(y=50,x=Year, label=N), vjust=1.2) +
    scale_y_continuous(expression(Coral~cover~('%'))) +
    scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
    theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
          plot.margin = margin(t=2,r=7,b=0,l=0),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.text.x=element_text(size=rel(1.2)),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color='gray70'),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
          strip.text=element_text(margin=margin(t=1, b=1, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
ggsave(file='output/figures/mantaCAT_centralNoBanner.png', g.central + facet_null(), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_centralNoBanner.pdf', g.central + facet_null(), width=5, height=3.5, units='in',dpi=300)

gt3=gt+geom_polygon(data=fortify(whagbr.c), aes(y=lat, x=long),fill=hues[4],color=NA)+
    geom_polygon(data=fortify(whagbr.c), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])
g <- ggplot_gtable(ggplot_build(g.central))
facets <- grep("strip-t-1-1", g$layout$name)
gg.central <- with(g$layout[facets,],
                   gtable_add_grob(g, ggplotGrob(gt3),t=t, l=5, b=b, r=6, name="pic_predator"))

ggsave(file='output/figures/mantaCAT_central.png', grid.draw(gg.central), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_central.pdf', grid.draw(gg.central), width=5, height=3.5, units='in',dpi=300)


## Southern GBR
g.southern<-newdata %>% filter(Location=='Southern GBR') %>% ggplot(aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2])+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd %>% filter(Location=='Southern GBR'), aes(y=50,x=Year, label=N), vjust=1.2) +
    scale_y_continuous(expression(Coral~cover~('%'))) +
    scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
    theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
          plot.margin = margin(t=2,r=7,b=0,l=0),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.text.x=element_text(size=rel(1.2)),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color='gray70'),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
          strip.text=element_text(margin=margin(t=1, b=1, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
ggsave(file='output/figures/mantaCAT_southernNoBanner.png', g.southern + facet_null(), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_southernNoBanner.pdf', g.southern + facet_null(), width=5, height=3.5, units='in',dpi=300)

gt4=gt+geom_polygon(data=fortify(whagbr.s), aes(y=lat, x=long),fill=hues[4],color=NA)+
    geom_polygon(data=fortify(whagbr.s), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])
g <- ggplot_gtable(ggplot_build(g.southern))
facets <- grep("strip-t-1-1", g$layout$name)
gg.southern <- with(g$layout[facets,],
                    gtable_add_grob(g, ggplotGrob(gt4),t=t, l=5, b=b, r=6, name="pic_predator"))

ggsave(file='output/figures/mantaCAT_southern.png', grid.draw(gg.southern), width=5, height=3.5, units='in',dpi=300)
ggsave(file='output/figures/mantaCAT_southern.pdf', grid.draw(gg.southern), width=5, height=3.5, units='in',dpi=300)



##Now generate the table==============================================================
load(file='data/modelled/dat.all.gbr.RData')
load(file='data/modelled/dat.all.northern.RData')
load(file='data/modelled/dat.all.central.RData')
load(file='data/modelled/dat.all.southern.RData')
newdata = rbind(dat.gbr,dat.northern,dat.central,dat.southern)
newdata$Location <- factor(newdata$Location, levels=unique(newdata$Location), labels=c('Great Barrier Reef','Northern GBR', 'Central GBR', 'Southern GBR'))

annualTable = newdata %>% dplyr:::select(-X1) %>% arrange(Location,Year)

annualTable =annualTable %>% filter(as.numeric(as.character(Year))>2010) %>% dplyr:::select(Location, Year, mean,lower, upper) %>% mutate_at(vars(mean,lower,upper), funs(round(.*100,1)))
dat.all = rbind(dat.all.gbr %>% ungroup %>% mutate(Location='Great Barrier Reef'),
                dat.all.northern,dat.all.central,dat.all.southern)
dat.all$Location <- factor(dat.all$Location, levels=unique(dat.all$Location), labels=c('Great Barrier Reef','Northern GBR', 'Central GBR', 'Southern GBR'))

a=dat.all %>% group_by(Location,Year) %>% summarize(N=n()) %>% mutate(Year=factor(Year)) %>% data.frame 
annualTable = annualTable %>% left_join(a)    


save(annualTable, file='output/tables/annualTable.RData')

# Now the full annual Table
fullannualTable = newdata %>% dplyr:::select(-X1) %>% arrange(Location,Year)

fullannualTable =fullannualTable %>% dplyr:::select(Location, Year, mean,lower, upper) %>% mutate_at(vars(mean,lower,upper), funs(round(.*100,1)))
dat.all = rbind(dat.all.gbr %>% ungroup %>% mutate(Location='Great Barrier Reef'),
                dat.all.northern,dat.all.central,dat.all.southern)
dat.all$Location <- factor(dat.all$Location, levels=unique(dat.all$Location), labels=c('Great Barrier Reef','Northern GBR', 'Central GBR', 'Southern GBR'))

a=dat.all %>% group_by(Location,Year) %>% summarize(N=n()) %>% mutate(Year=factor(Year)) %>% data.frame 
fullannualTable = fullannualTable %>% left_join(a)    


save(fullannualTable, file='output/tables/fullannualTable.RData')
