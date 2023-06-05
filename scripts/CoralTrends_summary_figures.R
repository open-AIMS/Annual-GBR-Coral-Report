source('CoralTrends_functions.R')
CoralTrends_checkPackages()
source('CoralTrends_config.R')

## ---- loadData
## original model
## load(file='../data/modelled/dat.northern_glmer.RData')
## load(file='../data/modelled/dat.central_glmer.RData')
## load(file='../data/modelled/dat.southern_glmer.RData')

## INLA reef, beta
## load(file='../data/modelled/mod.northern_inla_reef.beta.RData')
## load(file='../data/modelled/mod.central_inla_reef.beta.RData')
## load(file='../data/modelled/mod.southern_inla_reef.beta.RData')

## INLA reef, beta scaled
## load(file='../data/modelled/mod.northern_inla_reef.beta.scaled.RData')
## load(file='../data/modelled/mod.central_inla_reef.beta.scaled.RData')
## load(file='../data/modelled/mod.southern_inla_reef.beta.scaled.RData')

## glmmTMB reef, beta
## load(file='../data/modelled/dat.northern_glmmTMB.reef.beta.RData')
## load(file='../data/modelled/dat.central_glmmTMB.reef.beta.RData')
## load(file='../data/modelled/dat.southern_glmmTMB.reef.beta.RData')

## glmmTMB reef, beta disp
## load(file='../data/modelled/dat.northern_glmmTMB.reef.beta.disp.RData') #
## load(file='../data/modelled/dat.central_glmmTMB.reef.beta.disp.RData')  #
## load(file='../data/modelled/dat.southern_glmmTMB.reef.beta.disp.RData') #

## INLA tow, beta
## load(file='../data/modelled/mod.northern_inla.beta.RData') #
## load(file='../data/modelled/mod.central_inla.beta.RData')  #
## load(file='../data/modelled/mod.southern_inla.beta.RData') #

## INLA tow, beta disp
load(file='../data/modelled/mod.northern_inla.beta.disp.RData') #
dat.northern_inla.beta.disp <- dat.northern_inla.beta.disp %>%
    rename(response = mean)
load(file='../data/modelled/mod.central_inla.beta.disp.RData')  #
dat.central_inla.beta.disp <- dat.central_inla.beta.disp %>%
    rename(response = mean)
load(file='../data/modelled/mod.southern_inla.beta.disp.RData') #
dat.southern_inla.beta.disp <- dat.southern_inla.beta.disp %>%
    rename(response = mean)

## INLA tow, beta ry disp
load(file='../data/modelled/mod.northern_inla.beta.ry.disp.RData') #
dat.northern_inla.beta.ry.disp <- dat.northern_inla.beta.ry.disp %>%
    rename(response = mean)
load(file='../data/modelled/mod.central_inla.beta.ry.disp.RData')  #
dat.central_inla.beta.ry.disp <- dat.central_inla.beta.ry.disp %>%
    rename(response = mean)
load(file='../data/modelled/mod.southern_inla.beta.ry.disp.RData') #
dat.southern_inla.beta.ry.disp <- dat.southern_inla.beta.ry.disp %>%
    rename(response = mean)


## glmmTMB tow, beta
## load(file='../data/modelled/mod.northern_glmmTMB.beta.RData') #
## load(file='../data/modelled/mod.central_glmmTMB.beta.RData')  #
## load(file='../data/modelled/mod.southern_glmmTMB.beta.RData') #

## glmmTMB tow, beta disp
## load(file='../data/modelled/mod.gbr_glmmTMB.beta.disp.RData') #
load(file='../data/modelled/mod.northern_glmmTMB.beta.disp.RData') #
load(file='../data/modelled/mod.central_glmmTMB.beta.disp.RData')  #
load(file='../data/modelled/mod.southern_glmmTMB.beta.disp.RData') #

## glmmTMB tow, beta ry disp
## load(file='../data/modelled/mod.gbr_glmmTMB.beta.ry.disp.RData') #
load(file='../data/modelled/mod.northern_glmmTMB.beta.ry.disp.RData') #
load(file='../data/modelled/mod.central_glmmTMB.beta.ry.disp.RData')  #
load(file='../data/modelled/mod.southern_glmmTMB.beta.ry.disp.RData') #

## BRMS tow, beta disp
load(file='../data/modelled/mod.northern_brms.beta.disp.RData')
load(file='../data/modelled/mod.central_brms.beta.disp.RData')
load(file='../data/modelled/mod.southern_brms.beta.disp.RData')

## BRMS tow, beta ry disp
load(file='../data/modelled/mod.northern_brms.beta.ry.disp.RData')
load(file='../data/modelled/mod.central_brms.beta.ry.disp.RData')
load(file='../data/modelled/mod.southern_brms.beta.ry.disp.RData')

## BRMS tow, Ordinal
## load(file='../data/modelled/mod.northern_brms.cumulative.RData')
## load(file='../data/modelled/mod.central_brms.cumulative.RData')
## load(file='../data/modelled/mod.southern_brms.cumulative.RData')

## glmmTMB tow, beta disp random slopes
## load(file='../data/modelled/mod.northern_glmmTMB.beta.disp.rs.RData')
## load(file='../data/modelled/mod.central_glmmTMB.beta.disp.rs.RData')
## load(file='../data/modelled/mod.southern_glmmTMB.beta.disp.rs.RData')

load('../data/processed/manta.tow.RData')
load('../data/processed/manta.sum.RData')
load(file='../data/spatial/spatial_3Zone.RData')
load(file='../data/spatial/whagbr.RData')
load(file='../data/spatial/whagbr.n.RData')
load(file='../data/spatial/whagbr.c.RData')
load(file='../data/spatial/whagbr.s.RData')
## ----end

## ---- definitions
hues <- RColorBrewer::brewer.pal(4, "Blues")
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

## ----end

## ---- rawMeans
{
    ## ---- functions
    bfun <- function(x) {
        if (length(x$Cover)>3 & !all(x$Cover==x$Cover[1])) {
            betareg::betareg(Cover ~ 1, data=x) %>% coef() %>% `[[`(1) %>% plogis()
        } else {
            x$Cover[1]
        }
    }
    ## ----end

    ## ---- From Reef level data
    manta.stats.reef <- manta.sum %>%
        mutate(Year=REPORT_YEAR) %>% 
        group_by(Region, Year) %>%
        nest() %>%
        mutate(Mean=map_dbl(data, ~mean(.x$Cover)),
               Median=map_dbl(data, ~median(.x$Cover)),
               Logis=map_dbl(data, ~.x$Cover %>% gtools::logit() %>% mean() %>% plogis()),
               Beta=map_dbl(data, ~bfun(.x)) 
               ) %>%
        dplyr::select(-data) %>%
        unnest(cols=c()) %>%
        relocate(Mean, Beta,.after=last_col()) %>%
        arrange(Region, Year) %>%
        ungroup() %>%
        as.data.frame()
    ## ----end
    ## ---- From Reef level data
    manta.stats.tow <- manta.tow %>%
        group_by(Region, REEF_NAME, Year) %>%
        nest() %>%
        mutate(Mean=map_dbl(data, ~mean(.x$Cover)),
               Median=map_dbl(data, ~median(.x$Cover)),
               Logis=map_dbl(data, ~.x$Cover %>% gtools::logit() %>% mean() %>% plogis()),
               Beta=map_dbl(data, ~bfun(.x)) 
               ) %>%
        dplyr::select(-data) %>%
        unnest(cols=c()) %>%
        relocate(Mean, Beta,.after=last_col()) %>%
        arrange(Region, REEF_NAME, Year) %>%
        ungroup()

    manta.stats.tow.mean <- manta.stats.tow %>%
        mutate(Cover=Mean) %>% 
        group_by(Region, Year) %>%
        nest() %>%
        mutate(Mean=map_dbl(data, ~mean(.x$Cover)),
               Median=map_dbl(data, ~median(.x$Cover)),
               Logis=map_dbl(data, ~.x$Cover %>% gtools::logit() %>% mean() %>% plogis()),
               Beta=map_dbl(data, ~bfun(.x)) 
               ) %>%
        dplyr::select(-data) %>%
        unnest(cols=c()) %>%
        relocate(Mean, Beta,.after=last_col()) %>%
        arrange(Region, Year) %>%
        ungroup() %>%
        as.data.frame()

    manta.stats.tow.beta <- manta.stats.tow %>%
        mutate(Cover=Beta) %>% 
        group_by(Region, Year) %>%
        nest() %>%
        mutate(Mean=map_dbl(data, ~mean(.x$Cover)),
               Median=map_dbl(data, ~median(.x$Cover)),
               Logis=map_dbl(data, ~.x$Cover %>% gtools::logit() %>% mean() %>% plogis()),
               Beta=map_dbl(data, ~bfun(.x)) 
               ) %>%
        dplyr::select(-data) %>%
        unnest(cols=c()) %>%
        relocate(Mean, Beta,.after=last_col()) %>%
        arrange(Region, Year) %>%
        ungroup() %>%
        as.data.frame()
    ## ----end
}
## ----end


## Generate the banners
## ---- generateBanner

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
gt


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

## one with different spacing
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

## ----end
## ---- NumberOfReefs
nd <- manta.tow %>% group_by(Region) %>%
    summarise(Year=mean(range(as.numeric(as.character(Year)))),
              N=paste0('N=',length(unique(REEF_NAME)))) %>%
    bind_rows(manta.tow %>% summarise(Year=mean(range(as.numeric(as.character(Year)))),
                                      N=paste0('N=', length(unique(REEF_NAME)))) %>%
              mutate(Region='GBR'))
## ----end

## model_source = 'inla_reef.beta.scaled'      # INLA, Reef level, Beta scaled
## model_source = 'glmmTMB.reef.beta'          # glmmTMB, Reef level, Beta
## model_source = 'glmmTMB.reef.beta.disp'     # glmmTMB, Reef level, Beta disp

## model_source = 'inla.beta'                  # INLA, Tow level, Beta
## model_source = 'glmmTMB.beta'               # glmmTMB, Tow level, Beta
## model_source = 'glmmTMB.beta.disp'          # glmmTMB, Tow level, Beta disp
## model_source = 'inla.beta.disp'                # INLA, Tow level, Beta disp

## model_source = 'brms.beta.disp'             # BRMS, Tow level, Beta disp
model_source = 'brms.beta.ry.disp'             # BRMS, Tow level, Beta disp
## model_source = 'brms.cumulative'            # BRMS, Tow level, Ordinal

## model_source = 'glmmTMB.beta.disp.rs'
include_n=FALSE
include_gbr <- FALSE
## ---- threePanel
{
    ## ---- defineData
    ##dat.gbr <- sym(paste0('dat.gbr_',model_source))
    dat.northern <- sym(paste0('dat.northern_',model_source))
    dat.central <- sym(paste0('dat.central_',model_source))
    dat.southern <- sym(paste0('dat.southern_',model_source))

    newdata =
        #dat.gbr %>% eval %>% mutate(Region='GBR') %>%
        #rbind(dat.northern %>% eval %>% mutate(Region='Northern GBR')) %>%
        dat.northern %>% eval %>% mutate(Region='Northern GBR') %>%
        rbind(dat.central %>% eval %>% mutate(Region='Central GBR')) %>%
        rbind(dat.southern %>% eval %>% mutate(Region='Southern GBR')) %>%
        mutate(Region=factor(Region, levels=unique(Region))) %>%
        rename_with(recode, lower.HPD = 'lower', upper.HPD='upper',
                    lower.CL = 'lower', upper.CL = 'upper',
                    conf.low = 'lower', conf.high = 'upper',
                    mean='response', estimate='response') 
    if (!include_gbr) newdata <- newdata %>% filter(Region!='GBR') %>% droplevels()
    write_csv(newdata, file=paste0('../data/modelled/modelled_',model_source,'.csv'))
    ## ----end
    ## ---- PlotWithRibbons
    {
        ## ---- initalPlot
        g1<-ggplot(newdata, aes(y = response, x = as.numeric(as.character(Year))))+
            geom_blank(aes(y=0.10,x=1995))+geom_blank(aes(y=0.35,x=1995))+
            facet_wrap(~Region, nrow = 1, scales='fixed',
                       labeller = labeller(Region = setNames(paste0("\n", levels(newdata$Region),"\n"), levels(newdata$Region))))+
            geom_blank()+
            geom_ribbon(aes(ymin=lower, ymax=upper),fill=hues[2])+
                                        #geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
            geom_line(aes(x = as.numeric(as.character(Year))), color='blue') +
                                        #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
            scale_y_continuous(expression(Coral~cover~('%')),labels=function(x) x*100, expand=c(0,0),limits=c(0,0.50)) +
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
                  plot.margin=unit(c(0,0,2,0),'pt'),
                  panel.spacing.x=unit(10,'pt'))
        if (include_n)
            g1 <- g1 + geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) 
        g1
        ## ----end
        ## ---- addBanner
        if(!include_gbr) gt1=gt2; gt2=gt3; gt3=gt4;  
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
        grid.draw(gg)
        ## ----end
        ## ---- savePlot
        save(gg, file=paste0('../data/spatial/threePanels_',model_source,'.RData'))
        ggsave(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'), gg, width=15, height=3, units='in',dpi=300) 
        ggsave(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), gg, width=15, height=3, units='in',dpi=300) 
        png(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), width=15, height=3, units='in', res=300)
        grid.draw(gg)
        dev.off()
        ## ----end
        ## ---- addWatermark
        library(magick)
        library(png)
        a=magick::image_read(path='../parameters/ECM_1280945_v1_AIMS Logo Stacked white 1200px (2).png') #%>%
        gt1a=gt2 + annotation_custom(rasterGrob(a, x=unit(0.92,'npc'), y=unit(0.95, 'npc'), vjust=1,width=unit(0.1,'npc')))
        gt2a=gt3 + annotation_custom(rasterGrob(a, x=unit(0.92,'npc'), y=unit(0.95, 'npc'), vjust=1,width=unit(0.1,'npc')))
        gt3a=gt4 + annotation_custom(rasterGrob(a, x=unit(0.92,'npc'), y=unit(0.95, 'npc'), vjust=1,width=unit(0.1,'npc')))  
        gT <- ggplot_gtable(ggplot_build(g1))
        facets <- grep("strip-t-1-1", gT$layout$name)
        gg <- with(gT$layout[facets,],
                   gtable_add_grob(gT, ggplotGrob(gt1a),t=t, l=5, b=b, r=6, name="pic_predator"))
        facets <- grep("strip-t-2-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt2a),t=t, l=9, b=b, r=10, name="pic_predator"))
        facets <- grep("strip-t-3-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt3a),t=t, l=13, b=b, r=14, name="pic_predator"))
        grid.draw(gg)
        ggsave(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'), gg, width=15, height=3.5, units='in',dpi=300) 
        ggsave(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), gg, width=15, height=3.5, units='in',dpi=300) 
        png(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), width=15, height=3, units='in', res=300)
        grid.draw(gg)
        dev.off()

        ## In 2022, Comms insisted on a different logo for the
        ## watermark Unfortunately, it no longer fits in the strip.
        ## They also wanted it washed out
        a=magick::image_read(path='../parameters/AIMSLogo_Colour_inline.png') #%>%

        a <- a %>% magick::image_colorize(opacity = 50, color = "white") 
        g1b <- g1 + annotation_custom(rasterGrob(a, x=unit(0.05,'npc'), y=unit(0.3, 'npc'), vjust=1, hjust = 0, width=unit(0.4,'npc')))
        gT <- ggplot_gtable(ggplot_build(g1b))
        facets <- grep("strip-t-1-1", gT$layout$name)
        gg <- with(gT$layout[facets,],
                   gtable_add_grob(gT, ggplotGrob(gt2),t=t, l=5, b=b, r=6, name="pic_predator"))
        facets <- grep("strip-t-2-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=9, b=b, r=10, name="pic_predator"))
        facets <- grep("strip-t-3-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt4),t=t, l=13, b=b, r=14, name="pic_predator"))
        grid.draw(gg)
        ggsave(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'), gg, width=15, height=3.5, units='in',dpi=300) 
        ggsave(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), gg, width=15, height=3.5, units='in',dpi=300) 
        png(file=paste0('../output/figures/threePanels_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), width=15, height=3.5, units='in', res=300)
        grid.draw(gg)
        dev.off()
        ## ----end

    }
    ## ----end
    ## ---- PlotWithBars
    {
        ## ---- initalPlot
        g1<-ggplot(newdata, aes(y = response, x = as.numeric(as.character(Year))))+
            geom_blank(aes(y=0.10,x=1995))+geom_blank(aes(y=0.35,x=1995))+
            facet_wrap(~Region, nrow = 1, scales='fixed',
                       labeller = labeller(Region = setNames(paste0("\n", levels(newdata$Region),"\n"), levels(newdata$Region))))+
            geom_blank()+
            geom_pointrange(aes(ymin=lower, ymax=upper))+
            geom_line(aes(x = as.numeric(as.character(Year))), color='blue') +
            scale_y_continuous(expression(Coral~cover~('%')),labels=function(x) x*100, expand=c(0,0),limits=c(0,0.50)) +
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
                  strip.text=element_text(margin=margin(t=1, b=1,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.7,vjust=-1),
                  plot.margin=unit(c(0,0,2,0),'pt'),
                  panel.spacing.x=unit(10,'pt'))
        if (include_n)
            g1 <- g1 + geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) 
        g1
        ## ----end
        ## ---- addBanner
        if(!include_gbr) gt1=gt2; gt2=gt3; gt3=gt4;  
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
        grid.draw(gg)
        ## ----end
        ## ---- savePlot
        save(gg, file=paste0('../data/spatial/threePanels.Bars_',model_source,'.RData'))
        ggsave(file=paste0('../output/figures/threePanels.Bars_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'), gg, width=15, height=3, units='in',dpi=300) 
        ggsave(file=paste0('../output/figures/threePanels.Bars_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), gg, width=15, height=3, units='in',dpi=300) 
        png(file=paste0('../output/figures/threePanels.Bars_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), width=15, height=3, units='in', res=300)
        grid.draw(gg)
        dev.off()
        ## ----end
    }
    ## ----end
}
## ----end
## ---- singlePanels
same_y_axis_range <- TRUE
if (same_y_axis_range) {
    d <- sym(paste0('dat.northern_', model_source)) %>% eval() %>%
        bind_rows(sym(paste0('dat.central_', model_source)) %>% eval()) %>%
        bind_rows(sym(paste0('dat.southern_', model_source)) %>% eval()) %>%
        rename_with(recode, lower.HPD = 'lower', upper.HPD='upper',
                    lower.CL = 'lower', upper.CL = 'upper',
                    conf.low = 'lower', conf.high = 'upper',
                    mean='response', estimate='response') %>%
        ## rename(any_of(c(lower.CL = "lower", upper.CL = "upper"))) %>%
        summarise(y_range_min = min(lower),
                  y_range_max = max(upper))
}
## can also invoke region='GBR' - but only for glmmTMB.beta.disp
## for (region in c('GBR','Northern GBR', 'Central GBR', 'Southern GBR')) {
for (region in c('Northern GBR', 'Central GBR', 'Southern GBR')) {
    ## ---- prepareData
    nd1 <- nd %>% filter(Region==region)
    ## if (region=='GBR') nd1 <- manta.tow %>% summarise(Year=mean(range(as.numeric(as.character(Year)))),
    ##                                                   N=paste0('N=', length(unique(REEF_NAME))))
    reg <- dplyr::case_when(
                      ## region == 'GBR' ~ 'dat.gbr',
                      region == 'Northern GBR' ~ 'dat.northern',
                      region == 'Central GBR' ~'dat.central',
                      region == 'Southern GBR' ~ 'dat.southern')
    reg.shp <- dplyr::case_when(
                          ## region == 'GBR' ~ 'whagbr',
                          region == 'Northern GBR' ~ 'whagbr.n',
                          region == 'Central GBR' ~ 'whagbr.c',
                          region == 'Southern GBR' ~ 'whagbr.s'
                      ) %>%
        sym() %>%
        eval()
    dat <- sym(paste0(reg,'_',model_source)) %>% eval() %>%
        mutate(Region=region) %>%
        rename_with(recode, lower.HPD='lower', upper.HPD='upper',
                    lower.CL='lower', upper.CL='upper',
                    conf.low='lower', conf.high='upper',
                    mean='response', estimate='response')
    FILENAME=paste0('../data/modelled/modelled.gbr.',model_source)
    write_csv(dat, file=paste0(FILENAME,'.csv'))
    ## JS
    sink(paste0(FILENAME, '.js'))
    dat %>% send_df_to_js("data")
    sink()  

    ## ----end
    ## ---- PlotWithRibbons
    g <- dat %>%
        ggplot(aes(y=response, x=as.numeric(as.character(Year))))+
        geom_blank(aes(y=0.10,x=1995))+geom_blank(aes(y=0.35,x=1995))+
        facet_wrap(~Region, nrow=1, scales='fixed',
                   labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Region),"\n"), levels(newdata$Region))))+
        geom_blank()+
        geom_ribbon(aes(ymin=lower, ymax=upper),fill=hues[2])+
        geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
        scale_y_continuous(expression(Coral~cover~('%')), labels=function(x) x*100) +
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
              ## strip.text=element_text(margin=margin(t=1, b=1, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
              strip.text=element_text(margin=margin(t=1, b=1, r=10, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
    if (include_n)
        g <- g + geom_text(data=nd1, aes(y=0.50,x=Year, label=N), vjust=1.2)
    if (same_y_axis_range)
        g <- g + scale_y_continuous(expression(Coral~cover~('%')), labels=function(x) x*100, limits = c(0,0.5), expand = c(0,0)) 
    g

    # Watermarking
    library(magick)
    library(png)
    a=magick::image_read(path='../parameters/ECM_1280945_v1_AIMS Logo Stacked white 1200px (2).png') #%>%
    
    gt2=gt+geom_polygon(data=fortify(reg.shp), aes(y=lat, x=long),fill=hues[4],color=NA)+
        geom_polygon(data=fortify(reg.shp), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
        geom_polygon(fill='white', color=hues[4])
    ## Add watermark to the banner
    gt2 <- gt2 + annotation_custom(rasterGrob(a, x=unit(0.95,'npc'), y=unit(0.95, 'npc'), vjust=1,width=unit(0.1,'npc')))

    g <- ggplot_gtable(ggplot_build(g))
    facets <- grep("strip-t-1-1", g$layout$name)
    gg <- with(g$layout[facets,],
                        gtable_add_grob(g, ggplotGrob(gt2),t=t, l=5, b=b, r=6, name="pic_predator"))

    grid.draw(gg)

    ggsave(file=paste0('../output/figures/manta.',region,'_',model_source,'_',ifelse(include_n,'with_n',''),'.png'),
           gg, width=5, height=3.5, units='in',dpi=300)
    ggsave(file=paste0('../output/figures/manta.',region,'_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'),
           gg, width=5, height=3.5, units='in',dpi=300)
    ## ----end

    ## ---- PlotWithRibbons new logo
    g <- dat %>%
        ggplot(aes(y=response, x=as.numeric(as.character(Year))))+
        geom_blank(aes(y=0.10,x=1995))+geom_blank(aes(y=0.35,x=1995))+
        facet_wrap(~Region, nrow=1, scales='fixed',
                   labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Region),"\n"), levels(newdata$Region))))+
        geom_blank()+
        geom_ribbon(aes(ymin=lower, ymax=upper),fill=hues[2])+
        geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
        scale_y_continuous(expression(Coral~cover~('%')), labels=function(x) x*100) +
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
              ## strip.text=element_text(margin=margin(t=1, b=1, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
              strip.text=element_text(margin=margin(t=1, b=1, r=5, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
    if (include_n)
        g <- g + geom_text(data=nd1, aes(y=0.50,x=Year, label=N), vjust=1.2)
    if (same_y_axis_range)
        g <- g + scale_y_continuous(expression(Coral~cover~('%')), labels=function(x) x*100, limits = c(0,0.5), expand = c(0,0)) 
    g

    ## Now with the new logo
    library(magick)
    library(png)
    a=magick::image_read(path='../parameters/AIMSLogo_Colour_inline.png') #%>%
    ## a <- a %>% magick::image_colorize(opacity = 50, color = "white") 

    gt2=gt+geom_polygon(data=fortify(reg.shp), aes(y=lat, x=long),fill=hues[4],color=NA)+
        geom_polygon(data=fortify(reg.shp), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
        geom_polygon(fill='white', color=hues[4])
    
    gt2 <- gt2 + annotation_custom(rasterGrob(a, x=unit(0.85,'npc'), y=unit(0.9, 'npc'), vjust=1,width=unit(0.3,'npc')))
    ## g1b <- g1 + annotation_custom(rasterGrob(a, x=unit(0.05,'npc'), y=unit(0.3, 'npc'), vjust=1, hjust = 0, width=unit(0.4,'npc')))

    g <- ggplot_gtable(ggplot_build(g))
    facets <- grep("strip-t-1-1", g$layout$name)
    gg <- with(g$layout[facets,],
                        gtable_add_grob(g, ggplotGrob(gt2),t=t, l=5, b=b, r=6, name="pic_predator"))

    grid.draw(gg)

    ggsave(file=paste0('../output/figures/manta.',region,'_',model_source,'_',ifelse(include_n,'with_n',''),'.png'),
           gg, width=5, height=3.5, units='in',dpi=300)
    ggsave(file=paste0('../output/figures/manta.',region,'_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'),
           gg, width=5, height=3.5, units='in',dpi=300)

    ## ----end
    ## ---- PlotWithBars
    g <- dat %>%
        ggplot(aes(y=response, x=as.numeric(as.character(Year))))+
        geom_blank(aes(y=0.10,x=1995))+geom_blank(aes(y=0.35,x=1995))+
        facet_wrap(~Region, nrow=1, scales='fixed',
                   labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Region),"\n"), levels(newdata$Region))))+
        geom_blank()+
        geom_pointrange(aes(ymin=lower, ymax=upper))+
        geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
        scale_y_continuous(expression(Coral~cover~('%')), labels=function(x) x*100) +
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
              ## strip.text=element_text(margin=margin(t=1, b=1, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
              strip.text=element_text(margin=margin(t=1, b=1, r=5, unit='lines'),size=15,lineheight=0.5, face='bold',hjust=0.50,vjust=-1))
    if (include_n)
        g <- g + geom_text(data=nd1, aes(y=0.50,x=Year, label=N), vjust=1.2)
    g
    
    gt2=gt+geom_polygon(data=fortify(reg.shp), aes(y=lat, x=long),fill=hues[4],color=NA)+
        geom_polygon(data=fortify(reg.shp), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
        geom_polygon(fill='white', color=hues[4])
    ## Add watermark to the banner
    ## gt2 <- gt2 + annotation_custom(rasterGrob(a, x=unit(0.95,'npc'), y=unit(0.95, 'npc'), vjust=1,width=unit(0.1,'npc')))
    gt2 <- gt2 + annotation_custom(rasterGrob(a, x=unit(0.85,'npc'), y=unit(0.9, 'npc'), vjust=1,width=unit(0.3,'npc')))

    g <- ggplot_gtable(ggplot_build(g))
    facets <- grep("strip-t-1-1", g$layout$name)
    gg <- with(g$layout[facets,],
                        gtable_add_grob(g, ggplotGrob(gt2),t=t, l=5, b=b, r=6, name="pic_predator"))

    grid.draw(gg)
    
    ggsave(file=paste0('../output/figures/manta.',region,'_',model_source,'_',ifelse(include_n,'with_n',''),'_bars.png'),
           gg, width=5, height=3.5, units='in',dpi=300)
    ggsave(file=paste0('../output/figures/manta.',region,'_',model_source,'_',ifelse(include_n,'with_n',''),'_bars.pdf'),
           gg, width=5, height=3.5, units='in',dpi=300)
    ## ----end
}
## ----end

## ---- Zip
files <- list.files(path='../output/figures', pattern=paste0('^manta.*_',model_source), full.names=TRUE)
files1 <- list.files(path='../output/figures', pattern=paste0('^threePanels.*_',model_source,'(.p..|_with.*)'), full.names=TRUE)
files2 <- list.files(path='../data/modelled', pattern=paste0('^modelled_',model_source), full.names=TRUE)
files = c(files, files1, files2)
## files2 <- list.files(path='../output/figures', pattern=paste0('^threePanels.Bars.*_',model_source), full.names=TRUE)
zip(
    zipfile=paste0('../output/figures/Figures_',model_source,'_',ifelse(include_n,'with_n',''),'.zip'),
    files=files,
    flags='-rjo')

## ----end

## ---- Comparisons
RColorBrewer::brewer.pal(5, 'Set1')
mods <- tribble(
    ~model_source,              ~name,                             ~color,    ~linetype,
    'glmer',                    'Original (Reef level)',           '',        'solid',
    'inla_reef.beta',           'INLA Beta (Reef level)',          '',        'solid',
    'inla_reef.beta.scaled',    'INLA Beta scaled (Reef level)',   '',        'solid', 
    'glmmTMB.reef.beta',        'glmmTMB Beta (Reef level)',       '',        'solid',
    'glmmTMB.reef.beta.disp',   'glmmTMB Beta disp (Reef level)',  '',        'solid',
    'inla.beta',                'INLA Beta (Tow level)',           '',        'solid',
    'inla.beta.disp',           'INLA Beta disp (Tow level)',      '',        'solid',
    'inla.beta.ry.disp',        'INLA Beta ry disp (Tow level)',      '',        'solid',
    'glmmTMB.beta',             'glmmTMB Beta (Tow level)',        '',        'solid',
    'glmmTMB.beta.disp',        'glmmTMB Beta disp (Tow level)',   '#000000', 'solid',
    'glmmTMB.beta.ry.disp',        'glmmTMB Beta ry disp (Tow level)',   '#000000', 'solid',
    'brms.beta.disp',           'BRMS Beta disp (Tow level)',      '',        'solid',
    'brms.beta.ry.disp',           'BRMS Beta ry disp (Tow level)',      '',        'solid',
    'brms.cumulative',          'BRMS Cumulative (Tow level)',     '',        'solid',
    NA,                       'Raw Means',                         '',        'dashed',
    NA,                       'Raw Beta means',                    '',        'solid',
    )

newdata <- vector('list', nrow(mods %>% filter(!is.na(model_source))))
for (m in 1:nrow(mods %>% filter(!is.na(model_source)))) {
    print(mods$model_source[m])
    if (!exists(sym(paste0('dat.northern_', mods$model_source[m])))) next
    d.n <- sym(paste0('dat.northern_', mods$model_source[m]))
    d.c <- sym(paste0('dat.central_', mods$model_source[m]))
    d.s <- sym(paste0('dat.southern_', mods$model_source[m]))

    nd = d.n %>% eval %>% mutate(Region='Northern GBR') %>%
        rbind(d.c %>% eval %>% mutate(Region='Central GBR')) %>%
        rbind(d.s %>% eval %>% mutate(Region='Southern GBR')) %>%
        mutate(Region=factor(Region, levels=unique(Region))) %>%
        rename_with(recode, lower.HPD = 'lower', upper.HPD='upper',
                    lower.CL = 'lower', upper.CL = 'upper',
                    conf.low = 'lower', conf.high = 'upper',
                    mean='response', estimate='response')
    newdata[[m]] <- nd %>% mutate(Name=mods$name[m]) %>% dplyr::select(Name,Region,Year,response,lower,upper)
}
newdata = do.call('rbind', newdata)

newdata = newdata %>%
    rbind( manta.stats.reef %>% mutate(Name='Raw Means', lower=NA, upper=NA) %>% dplyr::select(Name,Region,Year,response=Mean,lower,upper)) %>%
    rbind( manta.stats.tow.beta %>% mutate(Name='Raw Beta means', lower=NA, upper=NA) %>% dplyr::select(Name,Region,Year,response=Mean,lower,upper)) 

## ---- Beta disp
{
    whichModels <- c(
        #'Raw Means',
        ##'Original (Reef level)',
#### 'INLA Beta (Reef level)',
#### 'INLA Beta scaled (Reef level)',
        #'glmmTMB Beta (Reef level)',
        ##'glmmTMB Beta disp (Reef level)',
        'INLA Beta disp (Tow level)',
        ## 'INLA Beta ry disp (Tow level)',
        ## 'glmmTMB Beta (Tow level)',
        'glmmTMB Beta disp (Tow level)',
        ## 'glmmTMB Beta ry disp (Tow level)',
        'BRMS Beta disp (Tow level)'#,
        #'BRMS Cumulative (Tow level)'
        ## 'Raw Beta means'
    )

    mods.used = mods %>% filter(name %in% whichModels) %>%
        mutate(N=1:n(),
               color=ifelse(color=='#000000', '#000000',
                            RColorBrewer::brewer.pal(max(N), 'Set1')[N]),
               fill=ifelse(color=='#000000','#000000',NA)
               ) %>%
        mutate(name=factor(name, levels=whichModels)) %>%
        arrange(name)

    g1 <- newdata %>%
        filter(Name %in% whichModels) %>%
        mutate(Name = factor(Name, levels=whichModels)) %>%
        ggplot() +
        geom_line(aes(y = response, x = as.numeric(as.character(Year)), colour = Name), alpha=1) +
        geom_ribbon(aes(y = response, ymin = lower, ymax = upper,
                        x = as.numeric(as.character(Year)),
                        fill = Name),
                    alpha=0.3) +
        facet_wrap(~Region, nrow=1) +
        ## scale_color_manual('Models', values=mods.used %>% pull(color)) +
        ## scale_linetype_manual('Models', values=mods.used %>% pull(linetype)) +
        ## scale_fill_manual('Models', values=mods.used %>% pull(fill)) +
        scale_x_continuous('') +
        scale_y_continuous('Hard coral cover (%)', labels=function(x) x*100) +
        theme_bw() +
        ## theme(legend.position='bottom',
        ##       legend.background = element_rect(fill='#00000020'),
        ##       legend.key=element_rect(fill='#00000000')) +
        theme(legend.position='bottom')+
        guides(colour = guide_legend(title = 'Model'),
               fill = guide_legend(title = 'Model'))
    ## guides(colour = guide_legend(title.position = "top", override.aes = list(fill=c(NA,NA,NA,NA,'#000000',NA,NA))))
    g1

    ggsave(file='../output/figures/ComparisonFigure_beta.disp.pdf', g1, width=9, height=3.5)
    ggsave(file='../output/figures/ComparisonFigure_beta.disp.png', g1, width=9, height=3.5, dpi=300)
}
## ----end
## ---- Beta ry disp
{
    whichModels <- c(
        #'Raw Means',
        ##'Original (Reef level)',
#### 'INLA Beta (Reef level)',
#### 'INLA Beta scaled (Reef level)',
        #'glmmTMB Beta (Reef level)',
        ##'glmmTMB Beta disp (Reef level)',
        ## 'INLA Beta disp (Tow level)',
        'INLA Beta ry disp (Tow level)',
        ## 'glmmTMB Beta (Tow level)',
        ## 'glmmTMB Beta disp (Tow level)',
        'glmmTMB Beta ry disp (Tow level)',
        ## 'BRMS Beta disp (Tow level)'#,
        'BRMS Beta ry disp (Tow level)'#,
        #'BRMS Cumulative (Tow level)'
        ## 'Raw Beta means'
    )

    mods.used = mods %>% filter(name %in% whichModels) %>%
        mutate(N=1:n(),
               color=ifelse(color=='#000000', '#000000',
                            RColorBrewer::brewer.pal(max(N), 'Set1')[N]),
               fill=ifelse(color=='#000000','#000000',NA)
               ) %>%
        mutate(name=factor(name, levels=whichModels)) %>%
        arrange(name)

    g1 <- newdata %>%
        filter(Name %in% whichModels) %>%
        mutate(Name = factor(Name, levels=whichModels)) %>%
        ggplot() +
        geom_line(aes(y = response, x = as.numeric(as.character(Year)), colour = Name), alpha=1) +
        geom_ribbon(aes(y = response, ymin = lower, ymax = upper,
                        x = as.numeric(as.character(Year)),
                        fill = Name),
                    alpha=0.3) +
        facet_wrap(~Region, nrow=1) +
        ## scale_color_manual('Models', values=mods.used %>% pull(color)) +
        ## scale_linetype_manual('Models', values=mods.used %>% pull(linetype)) +
        ## scale_fill_manual('Models', values=mods.used %>% pull(fill)) +
        scale_x_continuous('') +
        scale_y_continuous('Hard coral cover (%)', labels=function(x) x*100) +
        theme_bw() +
        ## theme(legend.position='bottom',
        ##       legend.background = element_rect(fill='#00000020'),
        ##       legend.key=element_rect(fill='#00000000')) +
        theme(legend.position='bottom')+
        guides(colour = guide_legend(title = 'Model'),
               fill = guide_legend(title = 'Model'))
    ## guides(colour = guide_legend(title.position = "top", override.aes = list(fill=c(NA,NA,NA,NA,'#000000',NA,NA))))
    g1

    ggsave(file='../output/figures/ComparisonFigure_beta.ry.disp.pdf', g1, width=9, height=3.5)
    ggsave(file='../output/figures/ComparisonFigure_beta.ry.disp.png', g1, width=9, height=3.5, dpi=300)
}
## ----end
## ---- Beta disp vs Beta ry disp
{
    whichModels <- c(
        #'Raw Means',
        ##'Original (Reef level)',
#### 'INLA Beta (Reef level)',
#### 'INLA Beta scaled (Reef level)',
        #'glmmTMB Beta (Reef level)',
        ##'glmmTMB Beta disp (Reef level)',
        'INLA Beta disp (Tow level)',
        'INLA Beta ry disp (Tow level)'#,
        ## 'glmmTMB Beta (Tow level)',
        ## 'glmmTMB Beta disp (Tow level)',
        ## 'glmmTMB Beta ry disp (Tow level)'#,
        ## 'BRMS Beta disp (Tow level)'#,
        ## 'BRMS Beta ry disp (Tow level)'#,
        #'BRMS Cumulative (Tow level)'
        ## 'Raw Beta means'
    )

    mods.used = mods %>% filter(name %in% whichModels) %>%
        mutate(N=1:n(),
               color=ifelse(color=='#000000', '#000000',
                            RColorBrewer::brewer.pal(max(N), 'Set1')[N]),
               fill=ifelse(color=='#000000','#000000',NA)
               ) %>%
        mutate(name=factor(name, levels=whichModels)) %>%
        arrange(name)

    g1 <- newdata %>%
        filter(Name %in% whichModels) %>%
        mutate(Name = factor(Name, levels=whichModels)) %>%
        ggplot() +
        geom_line(aes(y = response, x = as.numeric(as.character(Year)), colour = Name), alpha=1) +
        geom_ribbon(aes(y = response, ymin = lower, ymax = upper,
                        x = as.numeric(as.character(Year)),
                        fill = Name),
                    alpha=0.3) +
        facet_wrap(~Region, nrow=1) +
        ## scale_color_manual('Models', values=mods.used %>% pull(color)) +
        ## scale_linetype_manual('Models', values=mods.used %>% pull(linetype)) +
        ## scale_fill_manual('Models', values=mods.used %>% pull(fill)) +
        scale_x_continuous('') +
        scale_y_continuous('Hard coral cover (%)', labels=function(x) x*100) +
        theme_bw() +
        ## theme(legend.position='bottom',
        ##       legend.background = element_rect(fill='#00000020'),
        ##       legend.key=element_rect(fill='#00000000')) +
        theme(legend.position='bottom')+
        guides(colour = guide_legend(title = 'Model'),
               fill = guide_legend(title = 'Model'))
    ## guides(colour = guide_legend(title.position = "top", override.aes = list(fill=c(NA,NA,NA,NA,'#000000',NA,NA))))
    g1

    ggsave(file='../output/figures/ComparisonFigure_beta.disp.vs.beta.ry.disp.pdf', g1, width=9, height=3.5)
    ggsave(file='../output/figures/ComparisonFigure_beta.disp.vs.beta.ry.disp.png', g1, width=9, height=3.5, dpi=300)
}
## ----end


## ----end

## ---- left over rubbish
if(1==2) {

manta.tow %>% distinct(REEF_NAME, Year) %>%
    group_by(REEF_NAME) %>% count() %>%
    arrange(-n)

manta.stats <- manta.tow %>%
    ## filter(REEF_NAME== 'CHINAMAN REEF(22102)') %>%
    filter(REEF_NAME== '18023S') %>%
    group_by(Year) %>%
    nest() %>%
    mutate(Mean=map_dbl(data, ~mean(.x$Cover)),
           Median=map_dbl(data, ~median(.x$Cover)),
           Logis=map_dbl(data, ~.x$Cover %>% gtools::logit() %>% mean() %>% plogis()),
           Beta=map_dbl(data, ~ betareg::betareg(Cover ~ 1, data=.x) %>% coef() %>% `[[`(1) %>% plogis())
           ) %>%
    dplyr::select(-data) %>%
    unnest(cols=c()) %>%
    relocate(Mean, Beta,.after=last_col()) %>%
    arrange(Year) %>%
    as.data.frame()

manta.tow %>%
    filter(REEF_NAME=='CHINAMAN REEF(22102)') %>%
    ggplot() +
    geom_histogram(aes(x=Cover)) +
    geom_text(data=manta.stats, aes(y=Inf, x=Inf,
                                    label=paste0('Mean=',round(Mean*100,2),'\nBeta=',round(Beta*100,2)),
                                    hjust=1, vjust=1)) +
    scale_x_continuous('Cover', labels=function(x) x*100) +
    facet_wrap(~Year)

bfun <- function(x) {
    ## print(unique(x$RN))
    ## print(unique(x$YR))
    ## print(length(x$Cover))
    ## print(all(x$Cover==x$Cover[1]))
    ## print(x$Cover)
    if (length(x$Cover)>3 & !all(x$Cover==x$Cover[1])) {
        betareg::betareg(Cover ~ 1, data=x) %>% coef() %>% `[[`(1) %>% plogis()
    } else {
        x$Cover[1]
    }
}


manta.stats <- manta.tow %>%
    mutate(RN=REEF_NAME,YR=Year) %>%
    group_by(REEF_NAME, Year) %>%
    nest() %>%
    mutate(Mean=map_dbl(data, ~mean(.x$Cover)),
           Median=map_dbl(data, ~median(.x$Cover)),
           Logis=map_dbl(data, ~.x$Cover %>% gtools::logit() %>% mean() %>% plogis()),
           Beta=map_dbl(data, ~bfun(.x)) 
           ) %>%
    dplyr::select(-data) %>%
    unnest(cols=c()) %>%
    relocate(Mean, Beta,.after=last_col()) %>%
    arrange(Year) %>%
    ungroup() %>%
    as.data.frame()

manta.stats %>%
    group_by(Year) %>%
    summarise(Mean=mean(Mean),
              Beta=mean(Beta))

manta.stats %>%
    ggplot() +
    geom_abline(intercept=0, slope=1) +
    geom_point(aes(y=Beta, x=Mean))
}
## ----end
