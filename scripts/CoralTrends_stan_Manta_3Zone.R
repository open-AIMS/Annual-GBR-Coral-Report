source('CoralTrends_functions.R')
CoralTrends_checkPackages()
source('CoralTrends_config.R')
library(INLA)
library(emmeans)
library(DHARMa)
library(brms)

load('../data/processed/manta.sum.RData')
load('../data/processed/manta.tow.RData')

## Generate a list of reefs we are using to help accumulate other
## sources of associated data
all.reefs = manta.sum %>%
    dplyr:::select(P_CODE.mod,REEF_NAME,REEF_ID,Latitude,Longitude) %>%
    group_by(REEF_NAME,REEF_ID) %>%
    ## summarize_at(vars(Latitude,Longitude), funs(mean)) %>%
    summarize(across(c(Latitude,Longitude), mean)) %>%
    as.data.frame
write.csv(all.reefs,file='../data/all.reefs_3Zone.csv', quote=FALSE, row.names=FALSE)

## Genuine stan cannot handle proportional data for binomial families
## (particularly when weights are applied). A work-around is to
## multiple the proportion by the weights and convert this into an integer  
dat.all = manta.sum %>%
    mutate(Location=Region) %>%
    dplyr:::select(Cover, REEF_NAME, Tows,P_CODE.mod,Location,REPORT_YEAR) %>%
    mutate(Year=factor(REPORT_YEAR), N=length(unique(REEF_NAME))) %>% ungroup() %>%
    mutate(Cvr1 = as.integer(as.vector(Cover) * Tows), Cvr0 = Tows - Cvr1)


##original, stan_glmer beta, INLA_tow beta scaled, INLA_reef binomial, INLA_reef beta
##BRMS beta vanilla, BRMS beta disp, MGCV beta , MGCV ordinal, CLMM, BRMS_reef beta,
##
models <- c(
    ##'BRMS beta vanilla',
    #'BRMS beta disp'#,
    'BRMS beta ry disp'#,
    ##'glmmTMB beta vanilla',
    #'glmmTMB_tow beta disp',
    #'glmmTMB_tow beta ry disp'
    ##'BRMS ordinal'
)
zone <- c(
    ## 'northern'#,
    #'central'#,
    'southern'
)
COMPARE_MODELS <- FALSE
## Fit stan models===================================================================

if ( 1 == 2) {
## GBR
{
    ## ---- GBR
    {
        dat.all.gbr = dat.all %>% droplevels
        save(dat.all.gbr, file='../data/modelled/dat.all.gbr.RData')

        ## In 2021, we used the glmmTMB beta disp model


        ## ---- Gbr.Data
        {
            ## Tow level data
            manta.tow.gbr = manta.tow %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                                          levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                          ordered=TRUE),
                       nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL),
                       REEF_YEAR = interaction(REEF_NAME, Year)
                       )
            save(manta.tow.gbr, file='../data/modelled/manta.tow.gbr.RData')
        }
        ## ----end
        ## Raw cells ----------------------------------------------------------------
        ## ---- Raw cells
        {
            dat.all.gbr.cellmeans <- cellMeansRaw(dat.all.gbr)
            rawAdd <- ggproto_Raw(dat.all.gbr.cellmeans)

            manta.tow.gbr.cellmeans <- cellMeansRaw(manta.tow %>%
                                                    group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                                    summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                                    ungroup)
        }
        ## ----end                                        

        ## ---- GBR data prep
        {
            manta.tow.gbr = manta.tow %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                                          levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                          ordered=TRUE),
                       nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL),
                       REEF_YEAR = interaction(REEF_NAME, Year)
                       )
        }
        ## ----end

        ## ---- GBR.INLA.tow.ry.beta **
        {
            if ("INLA_tow beta ry disp" %in% models) {
                dat.inla <- dataINLA(dat=manta.tow.gbr, level='tow')
                dat = dat.inla[['dat.1']]
                dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
                rawAdd <- ggproto_Raw(dat.cellmeans)
                dat.scale = manta.tow.gbr %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                    summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
                ## we will leverage a mixed likelihood model
                dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
                
                dd <- dat %>%
                    dplyr::select(Cover, Year, REEF_NAME, REEF_YEAR) %>%
                    mutate(YEAR = Year)%>%
                    pivot_wider(id_cols = c(YEAR,REEF_NAME, REEF_YEAR),
                                names_from = Year,
                                values_from = Cover)
                dd1 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                    dplyr::select(matches("[0-9]{4}")) %>%
                    as.matrix()
                dd2 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                    dplyr::select(YEAR, REEF_NAME, REEF_YEAR)

                mod.gbr_inla.beta.ry.disp <- inla(form = dd1~YEAR +
                                                           f(REEF_NAME, model='iid') +
                                                           f(REEF_YEAR, model='iid'),
                                                       dat=dd2,
                                                       family=rep('beta',ncol(dd1)),
                                                       control.fixed = list(mean = 0, prec = 0.001,
                                                                            mean.intercept = 0.5,
                                                                            prec.intercept = 0.001),
                                                       control.predictor = list(compute = TRUE,
                                                                                link = 1,
                                                                                quantiles = c(0.025,0.25,0.5,0.75,0.975)
                                                                                )
                                                       )
                dat.gbr_inla.beta.ry.disp <- cellMeansINLA(mod=mod.gbr_inla.beta.ry.disp, newdata.hcc=dat.inla[['newdata.hcc']],
                                                                n.2=dat.inla[['n.2']], FUN=plogis)
                save(mod.gbr_inla.beta.ry.disp, dat.gbr_inla.beta.ry.disp, file='../data/modelled/mod.gbr_inla.beta.ry.disp.RData')
                save(dat.gbr_inla.beta.ry.disp, file='../data/modelled/dat.gbr_inla.beta.ry.disp.RData')
                rm(list=c('dat.gbr_inla.beta.ry.disp', 'mod.gbr_inla.beta.ry.disp'))
                gc()
            }
        }
        ## ----end
        ## ---- GBR.glmmTMB.tow.beta disp
        {
            if ('glmmTMB_tow beta disp' %in% models) {
                mod.gbr_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                                     dispformula = ~Year,
                                                     data=manta.tow.gbr,
                                                     ## weights=dat.all.gbr$Tows,
                                                     family=beta_family())
                dat.gbr_glmmTMB.beta.disp = emmeans(mod.gbr_glmmTMB.beta.disp, ~Year, type='response') %>%
                    as.data.frame()
                ## DHARMa::simulateResiduals(mod.gbr_glmmTMB.beta.disp, plot=TRUE)
                ## performance::check_model(mod.gbr_glmmTMB.beta.disp)
                save(mod.gbr_glmmTMB.beta.disp, dat.gbr_glmmTMB.beta.disp, file=paste0('../data/modelled/mod.gbr_glmmTMB.beta.disp.RData'))

            }
        }
        ## ----end
        ## ---- GBR.stan_glmer.reef.binomial
        {
            if ('original' %in% models) {
                mod.gbr <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                       data=dat.all.gbr, family=binomial,iter=5000,warmup=2500,
                                       chains=3,cores=3)
                dat.gbr <- data.frame(Location='Great Barrier Reef',Year=unique(dat.all.gbr$Year),
                                      N=length(unique(dat.all.gbr$REEF_NAME))) 
                Xmat <- model.matrix(~Year, dat.gbr)

                coefs = data.frame(mod.gbr) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
                Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
                dat.gbr = cbind(dat.gbr,
                                plyr:::adply(Fit,2,function(x) {
                                    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                                })
                                )
                save(dat.gbr, file='../data/modelled/dat.gbr.RData')
                save(mod.gbr, file='../data/modelled/mod.gbr.RData')

                l=levels(dat.all.gbr$Year)
                last_year=ranef(mod.gbr)[[1]][,c('(Intercept)',paste0('Year',l[(length(l)-1):length(l)]))]
                                        #last_year=ranef(mod.gbr)[[1]][,c('Year2016','Year2017')]
                ly = cbind(binomial()$linkinv(last_year[,1]-last_year[,2]),
                           binomial()$linkinv(last_year[,1]-last_year[,3])
                           )
                last_year$Diff = ly[,1]-ly[,2]
                last_year$REEF_NAME = rownames(last_year)
                save(last_year, file='../data/modelled/last_year.RData')
                
                rm(list=c('last_year','mod.gbr', 'dat.gbr', 'coefs', 'Fit'))
                gc()
            }
        }
        ## ----end
        ## ---- GBR.stan_glmer.reef.beta
        {
            if ("stan_glmer beta" %in% models) {
                mod.gbr <- stan_glmer(Cover ~ Year + (Year|P_CODE.mod/REEF_NAME),
                                      data=dat.all.gbr,
                                      family=mgcv::betar,
                                      iter=5000,
                                      warmup=2500,
                                      chains=3,cores=3,
                                      adapt_delta=0.95
                                      )
            }
        }
        ## ----end 
        ## ---- GBR.INLA.tow.beta
        {
            if ("INLA_tow beta scaled" %in% models) {
                dat.inla <- dataINLA(dat=manta.tow %>% droplevels,
                                     level='tow')
                dat = dat.inla[['dat.1']]
                dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
                rawAdd <- ggproto_Raw(dat.cellmeans)
                dat.scale = manta.tow %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                    summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
                mod.gbr <- ModelINLA_beta(form=Cover~Year +
                                              ## f(P_CODE.mod, model='iid'),#+
                                              ## f(REEF_NAME, model='iid') +
                                              ## f(REEF_NAME2, YEAR1, model='iid'),
                                              f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                              f(REEF_NAME2, YEAR1, model='iid', scale=dat.scale$Tows),
                                          ## f(REEF_NAME, model='iid', hyper=list(theta=list(initial=log(1), fixed=FALSE))) +
                                          ## f(REEF_NAME2, YEAR1, model='iid') ,
                                          ## f(REEF_NAME2, YEAR1, model='iid', hyper=list(theta=list(initial=log(1), fixed=FALSE, prior="pc.prec", param=c(3,0.05)))),
                                          ## f(REEF_NAME, YEAR, model='iid2d', n=2*length(unique(dat$REEF_NAME))) +
                                          ## f(REEF_NAME2, YEAR, copy="REEF_NAME"),
                                          dat=dat,
                                          family='beta',
                                          ## family='gaussian',
                                          weights=NULL)
                
                dat.gbr <- cellMeansINLA(mod=mod.gbr, newdata.hcc=dat.inla[['newdata.hcc']],
                                         n.2=dat.inla[['n.2']], FUN=plogis)
                dat.raw = dat %>% group_by(REEF_NAME, Year) %>% summarize(Cover=mean(Cover, na.rm=TRUE)) %>% filter(!is.na(Cover)) %>% droplevels
                dat.gbr %>%
                    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                    geom_line() +
                    rawAdd +
                    ggtitle('GBR INLA tow level')

                save(mod.gbr, file='../data/modelled/mod.gbr_inla_tow.RData')
                save(dat.gbr, file='../data/modelled/dat.gbr_inla_tow.RData')
                load(file='../data/modelled/mod.gbr_inla_tow.RData')
                load(file='../data/modelled/dat.gbr_inla_tow.RData')
                rm(list=c('dat.gbr', 'mod.gbr'))
                gc()
            }
        }
        ## ----end
        ## ---- GBR.INLA.tow.beta (REEF_YEAR)
        {
            if ("INLA_tow beta" %in% models) {
                dat.inla <- dataINLA(dat=manta.tow %>% droplevels,
                                     level='tow')
                dat = dat.inla[['dat.1']]
                dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
                rawAdd <- ggproto_Raw(dat.cellmeans)
                dat.scale = manta.tow %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                    summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
                mod.gbr <- ModelINLA_beta(form=Cover~Year +
                                              f(REEF_NAME, model='iid') +
                                              f(REEF_YEAR, model='iid'),
                                          dat=dat,
                                          family='beta',
                                          ## family='gaussian',
                                          weights=NULL)
                
                dat.gbr <- cellMeansINLA(mod=mod.gbr, newdata.hcc=dat.inla[['newdata.hcc']],
                                         n.2=dat.inla[['n.2']], FUN=plogis)
                dat.raw = dat %>% group_by(REEF_NAME, Year) %>% summarize(Cover=mean(Cover, na.rm=TRUE)) %>% filter(!is.na(Cover)) %>% droplevels
                dat.gbr %>%
                    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                    geom_line() +
                    rawAdd +
                    ggtitle('GBR INLA tow level')

                save(mod.gbr, file='../data/modelled/mod.gbr_inla_tow.RData')
                save(dat.gbr, file='../data/modelled/dat.gbr_inla_tow.RData')
                load(file='../data/modelled/mod.gbr_inla_tow.RData')
                load(file='../data/modelled/dat.gbr_inla_tow.RData')
                rm(list=c('dat.gbr', 'mod.gbr'))
                gc()
            }
        }
        ## ----end
        ## ---- GBR.INLA.reef.binomial
        {
            if ("INLA_reef binomial" %in% models) {
                dat.inla <- dataINLA(dat=dat.all.gbr %>% mutate(W=NA))
                dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
                mod.gbr <- ModelINLA_binomial(form=Cvr1~Year+
                                                  f(P_CODE.mod, model='iid')+
                                                  f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                                  f(REEF_NAME1, Year, model='iid', scale=dat.scale$Tows),
                                              dat=dat.inla[['dat.1']])
                newdata.gbr <- cellMeansINLA(mod=mod.gbr, newdata.hcc=dat.inla[['newdata.hcc']],
                                             n.2=dat.inla[['n.2']])
                newdata.gbr %>%
                    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                    geom_line() +
                    rawAdd
                save(mod.gbr, file='../data/modelled/mod.gbr_inla_binomial.RData')
                save(newdata.gbr, file='../data/modelled/newdata.gbr_inla_binomial.RData')
                rm(list=c('last_year','mod.gbr', 'newdata.gbr'))
                gc()
            }
        }
        ## ----end
        ## ---- GBR.INLA.reef.beta
        {
            if ("INLA_reef beta" %in% models) {
                dat.inla <- dataINLA(dat=dat.all.gbr %>% mutate(W=NA))
                dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
                mod.gbr <- ModelINLA_beta(form=Cover~Year+
                                              f(P_CODE.mod, model='iid')+
                                              f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                              f(REEF_NAME1, Year, model='iid', scale=dat.scale$Tows),
                                          dat=dat.inla[['dat.1']])
                newdata.gbr <- cellMeansINLA(mod=mod.gbr, newdata.hcc=dat.inla[['newdata.hcc']],
                                             n.2=dat.inla[['n.2']])
                newdata.gbr %>%
                    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                    geom_line() +
                    rawAdd
                save(mod.gbr,file='../data/modelled/mod.gbr_inla_beta.RData')
                save(newdata.gbr, file='../data/modelled/newdata.gbr_inla_beta.RData')
                rm(list=c('last_year','mod.gbr', 'newdata.gbr'))
                gc()
            }
        }
        ## ----end

        ## ---- Compare the models
        {
            load(file='../data/modelled/dat.gbr.RData')
            dat.gbr.original <- dat.gbr
            load(file='../data/modelled/mod.gbr.RData')
            load(file='../data/modelled/mod.gbr_inla_beta.RData')
            load(file='../data/modelled/newdata.gbr_inla_beta.RData')
            newdata.gbr_beta <- newdata.gbr
            load(file='../data/modelled/mod.gbr_inla_binomial.RData')
            load(file='../data/modelled/newdata.gbr_inla_binomial.RData')
            newdata.gbr_binomial <- newdata.gbr
            load(file='../data/modelled/dat.gbr_inla_tow.RData')

            ## original
            g1 <- dat.gbr.original %>%
                ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
                geom_line(color='blue') +
                scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
                scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
                scale_color_discrete('Raw data aggregate') + 
                rawAdd +
                theme_classic() +
                theme(axis.title.x=element_blank(),
                      legend.position=c(0.01,0.01), legend.justification=c(0,0),
                      panel.grid.minor=element_line(),
                      panel.grid.major=element_line()) +
                ggtitle('Original (stan binomial reef level)') +
                guides(color=guide_legend(nrow=2, byrow=TRUE))

            ## inla (reef level binomial)
            g2 <- newdata.gbr_binomial %>%
                ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
                geom_line(color='blue') +
                scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
                scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
                scale_color_discrete('Raw data aggregate') + 
                rawAdd +
                theme_classic() +
                theme(axis.title.x=element_blank(),
                      legend.position=c(0.01,0.01), legend.justification=c(0,0),
                      panel.grid.minor=element_line(),
                      panel.grid.major=element_line()) +
                ggtitle('INLA reef level binomial') +
                guides(color=guide_legend(nrow=2, byrow=TRUE))

            ## inla (reef level beta)
            g3 <- newdata.gbr_beta %>%
                ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
                geom_line(color='blue') +
                scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
                scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
                scale_color_discrete('Raw data aggregate') + 
                rawAdd +
                theme_classic() +
                theme(axis.title.x=element_blank(),
                      legend.position=c(0.01,0.01), legend.justification=c(0,0),
                      panel.grid.minor=element_line(),
                      panel.grid.major=element_line()) +
                ggtitle('INLA reef level beta') +
                guides(color=guide_legend(nrow=2, byrow=TRUE))

            ## inla (tow level beta)
            g4 <- dat.gbr %>%
                ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
                geom_line(color='blue') +
                scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
                scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
                scale_color_discrete('Raw data aggregate') + 
                rawAdd +
                theme_classic() +
                theme(axis.title.x=element_blank(),
                      legend.position=c(0.01,0.01), legend.justification=c(0,0),
                      panel.grid.minor=element_line(),
                      panel.grid.major=element_line()) +
                ggtitle('INLA tow level beta') +
                guides(color=guide_legend(nrow=2, byrow=TRUE))

            library(patchwork)
            g1 + g2 + g3 + g4


            rm(list=c('dat.gbr','mod.gbr','mod.gbr_inla_beta','newdata.gbr','newdata.gbr_beta', 'mod_gbr_inla_binomial','newdata.gbr_binomial'))

                                        #mod = glmmTMB(Cover ~ Year+(Year|P_CODE.mod/REEF_NAME), data=dat.all.gbr,
                                        #              family=beta_family)
        }
        ## ----end

        ## ---- Extract reef-level predictions (INLA.tow.beta)
        {
            mod.gbr$summary.fixed %>% head
            mod.gbr$summary.random[[1]] %>% dim
            mod.gbr$summary.random[[1]] %>% head
            reef_name <- manta.tow %>%  pull(REEF_NAME) %>% unique
            r1<-data.frame(REEF_NAME=reef_name, Rand.Int=mod.gbr$summary.random[[1]][,2])
            
            mod.gbr$summary.random[[2]] %>% dim
            mod.gbr$summary.random[[2]] %>% head
            ## reef_year <- manta.tow %>% droplevels %>% mutate(REEF_YEAR = paste0(REEF_NAME, Year)) %>% pull(REEF_YEAR) %>% unique
            ## reef_year %>% length
            reef_year <- dat[-dat.inla[['n.2']]] %>% pull(REEF_YEAR) %>% unique
            reef_year %>% length
            reef_year1=manta.tow %>% mutate(REEF_YEAR=factor(paste(REEF_NAME, Year))) %>% dplyr::select(REEF_NAME,Year,REEF_YEAR) %>% distinct %>% pull(REEF_YEAR)
            reef_year1 %>% length
            r2 <- data.frame(REEF_YEAR=reef_year[-length(reef_year)],
                             REEF_YEAR1=reef_year1[],
                             Rand.Slope=mod.gbr$summary.random[[2]][,2]) %>%
                mutate(REEF_NAME=gsub('(.*).[0-9]{4}','\\1',REEF_YEAR1))
            r3 <- r1 %>% left_join(r2)

            f1 <- mod.gbr$summary.fitted.values


            
            r3 %>% mutate(Fit=plogis(Rand.Int + Rand.Slope)) %>%
                head
        }
        ## ----end

        ## ---- Extract reef-level predictions (glmmTMB.tow.beta disp)
        {
            load(file=paste0('../data/modelled/mod.gbr_glmmTMB.beta.disp.RData'))

            ## coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_NAME['AGINCOURT REEFS (NO 1)',] 
            ## plogis(-1.068388 - 0.3291372)
            ## ndata = data.frame(Year=factor(1986:2021))
            ## Xmat <- model.matrix(~Year, ndata)
            ## coefs <- coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_NAME['AGINCOURT REEFS (NO 1)',] %>% as.matrix() %>% as.vector()
            ## ndata = cbind(ndata, Fit=as.vector(plogis(coefs %*% t(Xmat))))
            ## g1 + geom_line(data=ndata, aes(y=Fit, x=as.numeric(as.character(Year))), color='green')


            ## coefs <- coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_YEAR %>% filter(str_detect(rownames(.), 'AGINCOURT REEFS \\(NO 1\\)'))
            ## rownames(coefs)
            ## #pull(`(Intercept)`)
            ## coefs[1,]
            ## plogis(as.matrix(coefs[1,]) %*% t(Xmat))


            ## coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_YEAR['AGINCOURT REEFS (NO 1).1989:AGINCOURT REEFS (NO 1)',] 
            ## plogis(-1.514835 -0.3291372)




            ## plogis((-1.068388 + 0.3291372) + (-1.514835 +0.3291372))
            ## plogis((-1.068388 + 0.3291372) + (-1.514835 +0.3291372))


            ## aa =  fixef(mod.gbr_glmmTMB.beta.disp)[[1]]
            ## ndata = data.frame(Year=factor(1986:2021))
            ## Xmat <- model.matrix(~Year, ndata)

            ## ab1 = coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_YEAR %>%
            ##                                    rownames_to_column('R') %>%
            ##                                    filter(str_detect(R,'AGINCOURT REEFS \\(NO 1\\)'))
            ## ab1 %>% dim
            ## ab1 %>% head

            ## ab2 = coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_NAME %>%
            ##                                    rownames_to_column('R') %>%
            ##                                    filter(str_detect(R,'AGINCOURT REEFS \\(NO 1\\)'))
            ## ab2 %>% dim
            ## ab2 %>% head

            ## sweep(
            ##     x=as.matrix(ab1[,-1]),
            ##     MARGIN=2,
            ##     STATS=as.vector(as.matrix(ab2[-1])),
            ##     FUN='+'
            ##     )
            reefs <- manta.tow.gbr %>% pull(REEF_NAME) %>% unique
            fit.all.reefs <- vector('list', length(reefs))
            names(fit.all.reefs) <- reefs                     
            for (r in reefs) {
                print(r)
                raw.sum<-manta.tow.gbr %>% filter(REEF_NAME==r) %>%
                    group_by(Year) %>%
                    summarise(Cover=mean(Cover))

                                        #r='AGINCOURT REEFS (NO 1)'
                a0 <- fixef(mod.gbr_glmmTMB.beta.disp)[[1]][1]
                a1 <- fixef(mod.gbr_glmmTMB.beta.disp)[[1]][-1]
                a1 <- data.frame(Slope=a1) %>% rownames_to_column('Year') %>% mutate(Year=gsub('Year','',Year))
                a0 <- cbind(Intercept=a0, a1)

                a<-ranef(mod.gbr_glmmTMB.beta.disp) %>% `[[`(1) %>% `[[`(1) %>%
                    as.data.frame() %>%
                    rownames_to_column('R') %>%
                    mutate(REEF_NAME=gsub('(.*)\\.[0-9]{4}.*', '\\1', R),
                           Year=gsub('.*\\.([0-9]{4}).*', '\\1', R)) %>%
                    dplyr::select(-R) %>%
                    dplyr::rename(rand.Intercept=`(Intercept)`) %>%
                    filter(REEF_NAME==r)
                b<-ranef(mod.gbr_glmmTMB.beta.disp) %>% `[[`(1) %>% `[[`(2) %>%
                    as.data.frame %>%
                   rownames_to_column('REEF_NAME') %>%
                   dplyr::rename(rand.slope=`(Intercept)`) %>%
                   filter(REEF_NAME==r)
                a1 <- a %>% left_join(b)
                a3 <- a1 %>% left_join(a0) %>%
                    mutate(Cover=plogis(Intercept+rand.Intercept+Slope+rand.slope)) %>%
                    left_join(raw.sum %>% dplyr::rename(Raw=Cover))
                ## a3 %>% head
                fit.all.reefs[[r]] <- a3
                
                ## g1 <- ggplot() +
                ##     geom_line(data=a3, aes(y=Cover, x=as.numeric(as.character(Year))), color='red') +
                ##     geom_line(data=raw.sum, aes(y=Cover, x=as.numeric(as.character(Year))), color='blue')
                ## g1
            }
            fit.all.reefs <- do.call('rbind',fit.all.reefs)

            fit.all.reefs %>% filter(REEF_NAME=='BROOMFIELD REEF') %>%
                ggplot() +
                geom_line(aes(y=Cover, x=as.numeric(as.character(Year))),color='blue') +
                geom_line(aes(y=Raw, x=as.numeric(as.character(Year))), color='red')

            

            ## Using predict function is very slow and seems only to permit one prediction at a time with newdata?
            ## newdata <- data.frame(Year=1989:2021,
            ##                       REEF_NAME='AGINCOURT REEFS (NO 1)') %>%
            ##     mutate(REEF_YEAR=interaction(REEF_NAME,Year))
            
            ## predict(mod.gbr_glmmTMB.beta.disp, newdata=data.frame(Year=1990:1992, REEF_NAME=NA, REEF_YEAR=NA), re.form=~0, type='response')
            ## predict(mod.gbr_glmmTMB.beta.disp, newdata=newdata, re.form=NULL, type='response')
            ## predict(mod.gbr_glmmTMB.beta.disp, newdata=newdata[1,], type='response', se.fit=TRUE)



            ## predict(mod.gbr_glmmTMB.beta.disp, newdata=data.frame(Year=1990, REEF_NAME=r, REEF_YEAR=interaction(1990,r)))
        }
        ## ----end

        ## ---- GBR.glmmTMB.tow.beta.linear
        {
            mod.gbr_glmmTMB.beta.linear <- glmmTMB(Cover ~ REPORT_YEAR + (1|REEF_NAME/REEF_YEAR),
                                                   data=manta.tow.gbr,
                                                   ## weights=dat.all.gbr$Tows,
                                                   family=beta_family())
            dat.gbr_glmmTMB.beta.linear = emmeans(mod.gbr_glmmTMB.beta.linear, ~REPORT_YEAR, at=list(REPORT_YEAR=unique(manta.tow.gbr$REPORT_YEAR)), type='response') %>%
                as.data.frame()

            load(file=paste0('../data/modelled/mod.gbr_glmmTMB.beta.disp.RData'))
            dat.gbr_glmmTMB.beta.linear %>%
                ggplot() +
                geom_line(data=dat.gbr_glmmTMB.beta.disp, aes(y=response, x=as.numeric(as.character(Year))), color='blue') +
                geom_ribbon(data=dat.gbr_glmmTMB.beta.disp,aes(ymin=lower.CL, ymax=upper.CL, x=as.numeric(as.character(Year))), fill='lightblue', alpha=0.3) +
                geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, x=REPORT_YEAR), fill='orange', alpha=0.5) +
                geom_line(aes(y=response, x=REPORT_YEAR)) +
                scale_x_continuous('') +
                scale_y_continuous('Hard coral cover (%)', labels=function(x) x*100) +
                theme_classic() +
                ggtitle('GBR')

        }
        ## ----end
        
        ## ---- old
        {
            ##     ## Approximate the intercept per reef
            ##     r='AGINCOURT REEFS (NO 1)'
            ##     manta.tow.gbr %>% pull(REEF_NAME) %>% unique
            ##     raw.sum<-manta.tow.gbr %>% filter(reef_name==r) %>%
            ##         group_by(year) %>%
            ##         summarise(cover=mean(cover))
            ##     a0 <- fixef(mod.gbr_glmmTMB.beta.disp)[[1]][1]
            ##     a1 <- fixef(mod.gbr_glmmTMB.beta.disp)[[1]][-1]
            ##     a1 <- data.frame(Slope=a1) %>% rownames_to_column('Year') %>% mutate(Year=gsub('Year','',Year))
            
            ##     a<-ranef(mod.gbr_glmmTMB.beta.disp) %>% `[[`(1) %>% `[[`(1) %>%
            ##         as.data.frame %>%
            ##         ## mutate(Cover=plogis(`(Intercept)`)) %>%
            ##         rownames_to_column('R') %>%
            ##         mutate(REEF_NAME=gsub('(.*)\\.[0-9]{4}.*', '\\1', R),
            ##                Year=gsub('.*\\.([0-9]{4}).*', '\\1', R)) %>%
            ##         dplyr::select(-R) %>%
            ##         filter(REEF_NAME==r)
            ##     b<-ranef(mod.gbr_glmmTMB.beta.disp) %>% `[[`(1) %>% `[[`(2) %>%
            ##         as.data.frame %>%
            ##         rownames_to_column('REEF_NAME') %>%
            ##         mutate(Cover=`(Intercept)`) %>%
            ##         filter(REEF_NAME==r)
            ##         ## dplyr::select(-`(Intercept)`) %>%
            ##     a %>% left_join(b %>% dplyr::select(REEF_NAME, Cover)) %>%
            ##         mutate(C=a0) %>%
            ##         left_join(a1) %>% 
            ##         ## mutate(Cover1=plogis(C + `(Intercept)` + Cover)) %>%
            ##         mutate(Cover1=C + Cover,
            ##                Cover2 = `(Intercept)` + Slope,
            ##                Cover3 = plogis(Cover1 + Cover2)
            ##                ) %>%
            ##         ggplot() +
            ##         geom_line(aes(y=Cover3, x=as.numeric(as.character(Year))), color='red') +
            ##         geom_line(data=raw.sum, aes(y=Cover, x=as.numeric(as.character(Year))), color='blue')
            ## a %>% left_join(b %>% dplyr::select(REEF_NAME, Cover)) %>%
            ##     mutate(C=a0) %>%
            ##     left_join(a1) %>% 
            ##     ## mutate(Cover1=plogis(C + `(Intercept)` + Cover)) %>%
            ##     mutate(Cover1=C + Cover,
            ##            Cover2 = `(Intercept)` + Slope,
            ##            Cover3 = plogis(Cover1 + Cover2)
            ##            ) %>%
            ##     head   
        }
        ## ----end
    }
    ## ----end
}
}

## Northern
{
    ## ---- Northern
    ## ---- Northern.Data
    {
        dat.all.northern = dat.all %>%
            filter(Location=='Northern GBR') %>%
            droplevels %>%
            mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
            group_by(REEF_NAME) %>%
            mutate(W=mean(Tows, na.rm=TRUE)) %>%
            ungroup %>%
            mutate(W1=W/sum(W)) %>%
            group_by(Year) %>%
            mutate(W2=Tows/sum(Tows)) %>%
            ungroup
        save(dat.all.northern, file='../data/modelled/dat.all.northern.RData')
        ## Tow level data
        manta.tow.northern = manta.tow %>%
            filter(Region=='Northern GBR') %>%
            droplevels %>%
            mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                                      levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                      ordered=TRUE),
                   nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
                   nLIVE_CORAL=as.numeric(oLIVE_CORAL),
                   REEF_YEAR = interaction(REEF_NAME, Year)
                   )
        save(manta.tow.northern, file='../data/modelled/manta.tow.northern.RData')
    }
    ## ----end

    ## Raw cells
    ## ---- Northern.Raw cells
    {
        dat.all.northern.cellmeans <- cellMeansRaw(dat.all.northern)
        rawAdd <- ggproto_Raw(dat.all.northern.cellmeans)

        manta.tow.northern.cellmeans <- cellMeansRaw(manta.tow %>%
                                                     filter(Region=='Northern GBR') %>%
                                                     droplevels %>% 
                                                     group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                                     summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                                     ungroup)
    }
    ## ----end
    
    ## ---- Northern.stan_glmer.reef.binomial
    {
        if ('original' %in% models) {
            ## mod.northern <- ModelOriginal(form=formula(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME)),
            ##                               dat=dat.all.northern,
            ##                               location='Northern')
            ## mod.northern[['newdata']] %>%
            ##     ggplot(aes(y=mean, x=as.numeric(Year))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line()
            ## mod.northern[['mod']] %>% save(file='../data/modelled/mod.northern_original.RData')
            ## mod.northern[['newdata']] %>% save(file='../data/modelled/newdata.northern_original.RData')

            mod.northern_glmer <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                              data=dat.all.northern, family=binomial,iter=5000,warmup=2500,
                                              chains=3,cores=3)
            dat.northern_glmer <- data.frame(Location='Northern GBR',Year=unique(dat.all.northern$Year), N=length(unique(dat.all.northern$REEF_NAME)))
            Xmat <- model.matrix(~Year, dat.northern_glmer)

            coefs = data.frame(mod.northern_glmer) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
            Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
            dat.northern_glmer = cbind(dat.northern_glmer,
                                       plyr:::adply(Fit,2,function(x) {
                                           data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                                       })
                                       )
            save(dat.northern_glmer, mod.northern_glmer, file='../data/modelled/dat.northern_glmer.RData')
            ## save(dat.northern, file='../data/modelled/dat.northern.RData')
            ## save(mod.northern, file='../data/modelled/mod.northern.RData')
            rm(list=c('last_year','dat.northern_glmer','mod.northern_glmer'))
            gc()
            ## mod.northern <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
            ##                             data=dat.all.northern,
            ##                             family=binomial,
            ##                             iter=5000,
            ##                             warmup=2500,
            ##                             chains=3,cores=3)
        }
    }
    ## ----end
    ## ---- Northern.stan_glmer.reef.beta
    {
        if ('stan_glmer beta' %in% models) {
            mod.northern <- Modelstan_glmer(form=Cover ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                            dat=dat.all.northern)
            newdata.northern <- cellMeansOriginal(mod.northern, dat=dat.all.northern, location='Northern') 
            newdata.northern %>%
                ggplot(aes(y=mean, x=as.numeric(Year))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                geom_line()
            save(mod.northern, file='../data/modelled/mod.northern_stan_glmer.RData')
            save(newdata.northern, file='../data/modelled/newdata.northern_stan_glmer.RData')
            rm(list=c('last_year','mod.northern', 'newdata.northern', 'coefs', 'Fit'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.BRMS.reef.beta
    {
        if ('BRMS_reef beta' %in% models) {
            priors <- prior(normal(0, 5), class = "b") +
                prior(normal(0, 5), class = "Intercept") +
                ## prior(gamma(2, 0.1), class = "sd") +
                prior(gamma(1, 0.5), class = "sd") +
                prior(gamma(0.01, 0.01), class = "phi")
            inits = list(list(phi=list(rgamma(1,0.1,0.1))),
                         list(phi=list(rgamma(1,0.1,0.1))),
                         list(phi=list(rgamma(1,0.1,0.1)))
                         )
            mod.northern <- Modelbrms(form=Cover|weights(Tows) ~ Year + (1|REEF_NAME),
                                      dat=dat.all.northern)
            
            mod.northern <- Modelbrms(form=Cover|weights(Tows) ~ Year + (Year|REEF_NAME),
                                      dat=dat.all.northern)
        }
    }
    ## ----end
    ## ---- Northern.INLA.reef.beta vanilla
    {
        if ("INLA_reef beta" %in% models) {
            dat.inla <- dataINLA(dat=dat.all.northern %>% mutate(W=NA))
            dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
            mod.northern_inla_reef.beta <- ModelINLA_beta(form=Cover~Year+
                                                              f(P_CODE.mod, model='iid')+
                                                              f(REEF_NAME, model='iid'),
                                                          dat=dat.inla[['dat.1']])
            dat.northern_inla_reef.beta <- cellMeansINLA(mod=mod.northern_inla_reef.beta, newdata.hcc=dat.inla[['newdata.hcc']],
                                                         n.2=dat.inla[['n.2']])
            save(dat.northern_inla_reef.beta, mod.northern_inla_reef.beta, file='../data/modelled/mod.northern_inla_reef.beta.RData')
            rm(list=c('mod.northern_inla_reef.beta', 'dat.northern_inla_reef.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.INLA.reef.beta scaled
    {
        if ("INLA_reef beta" %in% models) {
            dat.inla <- dataINLA(dat=dat.all.northern %>% mutate(W=NA))
            dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
            mod.northern_inla_reef.beta.scaled <- ModelINLA_beta(form=Cover~Year+
                                                                     f(P_CODE.mod, model='iid')+
                                                                     f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                                                     f(REEF_NAME1, Year, model='iid', scale=dat.scale$Tows),
                                                                 dat=dat.inla[['dat.1']])
            dat.northern_inla_reef.beta.scaled <- cellMeansINLA(mod=mod.northern_inla_reef.beta.scaled, newdata.hcc=dat.inla[['newdata.hcc']],
                                                                n.2=dat.inla[['n.2']])
            ## newdata.northern %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd
            save(dat.northern_inla_reef.beta.scaled, mod.northern_inla_reef.beta.scaled, file='../data/modelled/mod.northern_inla_reef.beta.scaled.RData')
            ## save(newdata.northern, file='../data/modelled/newdata.northern_inla_beta.RData')
            rm(list=c('mod.northern_inla_reef.beta.scaled', 'dat.northern_inla_reef.beta.scaled'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.INLA.reef.binomial
    {
        if ('INLA_reef binomial' %in% models) {
            dat.inla <- dataINLA(dat=dat.all.northern)
            mod.northern <- ModelINLA_binomial(form=Cvr1~Year+
                                                   f(P_CODE.mod, model='iid')+
                                                   f(REEF_NAME, model='iid') +
                                                   f(REEF_NAME1, Year, model='iid'),
                                               dat=dat.inla[['dat.1']])
            newdata.northern <- cellMeansINLA(mod=mod.northern, newdata.hcc=dat.inla[['newdata.hcc']],
                                              n.2=dat.inla[['n.2']])
            newdata.northern %>%
                ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                geom_line() +
                rawAdd
            save(mod.northern, file='../data/modelled/mod.northern_inla_binomial.RData')
            save(newdata.northern,file='../data/modelled/newdata.northern_inla_binomial.RData')
            rm(list=c('mod.northern', 'newdata.northern'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.glmmTMB.reef.beta
    {
        if ('glmmTMB_reef beta' %in% models) {
            mod.northern_glmmTMB.reef.beta <- glmmTMB(Cover ~ Year + (1|P_CODE.mod) + (1|REEF_NAME),
                                                      data=dat.all.northern,
                                                      weights=dat.all.northern$Tows,
                                                      family=beta_family())
            dat.northern_glmmTMB.reef.beta  = emmeans(mod.northern_glmmTMB.reef.beta, ~Year, type='response') %>%
                as.data.frame()
            save(mod.northern_glmmTMB.reef.beta, dat.northern_glmmTMB.reef.beta, file='../data/modelled/dat.northern_glmmTMB.reef.beta.RData')
            
            ## emmeans(mod.northern, ~Year, type='response') %>%
            ##     as.data.frame() %>% 
            ##     ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd

            rm(list=c('dat.northern_glmmTMB.reef.beta', 'mod.northern_glmmTMB.reef.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.glmmTMB.reef.beta.disp
    {
        if ('glmmTMB_reef beta disp' %in% models) {
            mod.northern_glmmTMB.reef.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME),
                                                           dispformula = ~ Year,
                                                           data=dat.all.northern,
                                                           weights=dat.all.northern$Tows,
                                                           family=beta_family())
            dat.northern_glmmTMB.reef.beta.disp  = emmeans(mod.northern_glmmTMB.reef.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            save(mod.northern_glmmTMB.reef.beta.disp, dat.northern_glmmTMB.reef.beta.disp, file='../data/modelled/dat.northern_glmmTMB.reef.beta.disp.RData')
            
            ## emmeans(mod.northern, ~Year, type='response') %>%
            ##     as.data.frame() %>% 
            ##     ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd

            rm(list=c('dat.northern_glmmTMB.reef.beta.disp', 'mod.northern_glmmTMB.reef.beta.disp'))
            gc()
        }
    }
    ## ----end

    ## ---- Northern.INLA.tow.beta scaled
    {
        if ("INLA_tow beta scaled" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.northern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.northern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            mod.northern_inla.beta.scaled <- ModelINLA_beta(form=Cover~Year +
                                                                f(P_CODE.mod, model='iid') +
                                                                f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                                                f(REEF_NAME2, YEAR1, model='iid', scale=dat.scale$Tows),
                                                            dat=dat,
                                                            family='beta',
                                                            weights=NULL)
            
            dat.northern_inla.beta.scaled <- cellMeansINLA(mod=mod.northern_inla.beta.scaled, newdata.hcc=dat.inla[['newdata.hcc']],
                                                           n.2=dat.inla[['n.2']], FUN=plogis)
            ## dat.raw = dat %>% group_by(REEF_NAME, Year) %>% summarize(Cover=mean(Cover, na.rm=TRUE)) %>% filter(!is.na(Cover)) %>% droplevels
            ## dat.northern %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill='blue') +
            ##     geom_line(color='blue') +
            ##    rawAdd +
            ##    ggtitle('GBR INLA tow level') #
            ##    ## geom_line(data=mod.brms.df, aes(y=response, x=as.numeric(as.character(year)))) +
            ##    ## geom_ribbon(data=mod.brms.df, aes(y=response, ymin=lower.HPD, ymax=upper.HPD, x=as.numeric(as.character(year))), alpha=0.3)        

            save(mod.northern_inla.beta.scaled, dat.northern_inla.beta.scaled, file='../data/modelled/mod.northern_inla.beta.scaled.RData')
            rm(list=c('dat.northern_inla.beta.scaled', 'mod.northern_inla.beta.scaled'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.INLA.tow.beta
    {
        if ("INLA_tow beta" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.northern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.northern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            mod.northern_inla.beta <- ModelINLA_beta(form=Cover~Year +
                                                         f(P_CODE.mod, model='iid') +
                                                         f(REEF_NAME, model='iid') +
                                                         f(REEF_YEAR, model='iid'), 
                                                     dat=dat,
                                                     family='beta',
                                                     weights=NULL)
            dat.northern_inla.beta <- cellMeansINLA(mod=mod.northern_inla.beta, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.northern_inla.beta, dat.northern_inla.beta, file='../data/modelled/mod.northern_inla.beta.RData')
            rm(list=c('dat.northern_inla.beta', 'mod.northern_inla.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.INLA.tow.beta **
    {
        if ("INLA_tow beta disp" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.northern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.northern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            ## we will leverage a mixed likelihood model
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            
            dd <- dat %>%
                dplyr::select(Cover, Year, REEF_NAME) %>%
                mutate(YEAR = Year)%>%
                pivot_wider(id_cols = c(YEAR,REEF_NAME),
                            names_from = Year,
                            values_from = Cover)
            dd1 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(matches("[0-9]{4}")) %>%
                as.matrix()
            dd2 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(YEAR, REEF_NAME)

            mod.northern_inla.beta.disp <- inla(form = dd1~YEAR +
                                               f(REEF_NAME, model='iid'),
                                           dat=dd2,
                                           family=rep('beta',ncol(dd1)),
                                           control.fixed = list(mean = 0, prec = 0.001,
                                                                mean.intercept = 0.5,
                                                                prec.intercept = 0.001),
                                           control.predictor = list(compute = TRUE,
                                                                    link = 1,
                                                                    quantiles = c(0.025,0.25,0.5,0.75,0.975)
                                                                    )
                                           )
            dat.northern_inla.beta.disp <- cellMeansINLA(mod=mod.northern_inla.beta.disp, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.northern_inla.beta.disp, dat.northern_inla.beta.disp, file='../data/modelled/mod.northern_inla.beta.disp.RData')
            rm(list=c('dat.northern_inla.beta.disp', 'mod.northern_inla.beta.disp'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.INLA.tow.ry.beta **
    {
        if ("INLA_tow beta ry disp" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.northern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.northern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            ## we will leverage a mixed likelihood model
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            
            dd <- dat %>%
                dplyr::select(Cover, Year, REEF_NAME, REEF_YEAR) %>%
                mutate(YEAR = Year)%>%
                pivot_wider(id_cols = c(YEAR,REEF_NAME, REEF_YEAR),
                            names_from = Year,
                            values_from = Cover)
            dd1 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(matches("[0-9]{4}")) %>%
                as.matrix()
            dd2 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(YEAR, REEF_NAME, REEF_YEAR)

            mod.northern_inla.beta.ry.disp <- inla(form = dd1~YEAR +
                                               f(REEF_NAME, model='iid') +
                                               f(REEF_YEAR, model='iid'),
                                           dat=dd2,
                                           family=rep('beta',ncol(dd1)),
                                           control.fixed = list(mean = 0, prec = 0.001,
                                                                mean.intercept = 0.5,
                                                                prec.intercept = 0.001),
                                           control.predictor = list(compute = TRUE,
                                                                    link = 1,
                                                                    quantiles = c(0.025,0.25,0.5,0.75,0.975)
                                                                    )
                                           )
            dat.northern_inla.beta.ry.disp <- cellMeansINLA(mod=mod.northern_inla.beta.ry.disp, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.northern_inla.beta.ry.disp, dat.northern_inla.beta.ry.disp, file='../data/modelled/mod.northern_inla.beta.ry.disp.RData')
            rm(list=c('dat.northern_inla.beta.ry.disp', 'mod.northern_inla.beta.ry.disp'))
            gc()
        }
    }
    ## ----end
    ## ---- Northern.BRMS.tow.beta vanilla
    {
        if ('BRMS beta vanilla' %in% models) {
            mod.northern_brms.beta <- brm(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                          data=manta.tow.northern,
                                          family=Beta(link='logit'),
                                          iter=1e4,
                                          warmup=5e3,
                                          thin=5,
                                          chains=4, cores=4,
                                          prior = prior(normal(0, 3), class = "b") +
                                              prior(normal(0, 3), class = "Intercept") +
                                              prior(gamma(2, 1), class = "sd") +
                                              prior(gamma(2, 1), class = "sd", group = "REEF_NAME") +
                                              prior(gamma(2, 1), class = "phi")
                                          )
            dat.northern_brms.beta = emmeans(mod.northern_brms.beta, ~Year, type='response') %>%
                as.data.frame()
            save(mod.northern_brms.beta, dat.northern_brms.beta, file=paste0('../data/modelled/mod.northern_brms_beta.RData'))
            ## dat.northern_brms.beta %>%
            ##     ggplot(aes(y=response, x=as.numeric(as.character(year)))) +
            ##     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), alpha=0.3) +
            ##     geom_line() +
            ##     geom_ribbon(data=dat.northern, aes(y=mean, x=as.numeric(as.character(Year)), ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.northern, aes(y=mean, x=as.numeric(as.character(Year))), color='blue') +
            ##     geom_ribbon(data=dat.northern.original, aes(y=mean, x=as.numeric(as.character(Year)), ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.northern.original, aes(y=mean, x=as.numeric(as.character(Year))), color='red')
            ## save(dat.northern_brms.beta, file='../data/modelled/dat.northern_brms.beta.RData')
            rm(list=c('dat.northern_brms.beta', 'mod.northern_brms.beta'))
            gc()
        }    
    }
    ## ----end
    ## ---- Northern.BRMS.tow.beta disp **
    {
        if ('BRMS beta disp' %in% models & 'northern' %in% zone) {
            cat('Fitting brms disp for Northern\n\n')
            mod.northern_brms.beta.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME), phi~0+Year),
                                               data=manta.tow.northern,
                                               family=Beta(link='logit'),
                                               iter=1e4,
                                               warmup=5e3,
                                               thin=5,
                                               chains=4, cores=4,
                                               prior = prior(normal(0, 3), class = "b") +
                                                   prior(normal(0, 3), class = "Intercept") +
                                                   prior(gamma(2, 1), class = "sd") #+
                                               ## prior(gamma(2, 1), class = "phi")
                                               )
            ## ---- Northern.BRMS.tow.beta disp diagnostics
            {
                ## sampling diagnostics
                pdf(file = '../output/figures/traceplots_northern_brms.beta.disp.pdf')
                rstan::traceplot(mod.northern_brms.beta.disp$fit)
                dev.off()

                ## density overlay
                pdf(file = '../output/figures/density_northern_brms.beta.disp.pdf')
                mod.northern_brms.beta.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100) 
                dev.off()

                ## DHARMa residuals
                preds <- mod.northern_brms.beta.disp %>%
                    posterior_predict(nsamples = 250, summary = FALSE)
                mod.resids <- createDHARMa(
                    simulatedResponse = t(preds),
                    observedResponse = manta.tow.northern$Cover,
                    fittedPredictedResponse = apply(preds, 2, median),
                    integerResponse = FALSE
                )
                pdf(file = '../output/figures/DHARMa_northern_brms.beta.disp.pdf')
                mod.resids %>% plot()
                dev.off()
                save(mod.resids, file=paste0('../data/modelled/resids.northern_brms.beta.disp.RData'))
            }
            ## ----end
            dat.northern_brms.beta.disp = emmeans(mod.northern_brms.beta.disp,
                                                  ~Year, type='response') %>%
                as.data.frame()
            save(mod.northern_brms.beta.disp, dat.northern_brms.beta.disp, file=paste0('../data/modelled/mod.northern_brms.beta.disp.RData'))
            rm(list=c('dat.northern_brms.beta.disp', 'mod.northern_brms.beta.disp'))
            gc()
        }  
    }
    ## ----end
    ## ---- Northern.BRMS.tow.beta ry disp **
    {
        if ('BRMS beta ry disp' %in% models & 'northern' %in% zone) {
            cat('Fitting brms ry disp for Northern\n\n')
            mod.northern_brms.beta.ry.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME/REEF_YEAR), phi~0+Year),
                                               data=manta.tow.northern,
                                               family=Beta(link='logit'),
                                               iter=1e4,
                                               warmup=5e3,
                                               thin=5,
                                               chains=4, cores=4,
                                               prior = prior(normal(0, 3), class = "b") +
                                                   prior(normal(0, 3), class = "Intercept") +
                                                   prior(gamma(2, 1), class = "sd") #+
                                               ## prior(gamma(2, 1), class = "phi")
                                               )
            ## ---- Northern.BRMS.tow.beta disp diagnostics
            {
                ## sampling diagnostics
                pdf(file = '../output/figures/traceplots_northern_brms.beta.ry.disp.pdf')
                rstan::traceplot(mod.northern_brms.beta.ry.disp$fit)
                dev.off()

                ## density overlay
                pdf(file = '../output/figures/density_northern_brms.beta.ry.disp.pdf')
                mod.northern_brms.beta.ry.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100) 
                dev.off()

                ## DHARMa residuals
                preds <- mod.northern_brms.beta.ry.disp %>%
                    posterior_predict(nsamples = 250, summary = FALSE)
                mod.resids <- createDHARMa(
                    simulatedResponse = t(preds),
                    observedResponse = manta.tow.northern$Cover,
                    fittedPredictedResponse = apply(preds, 2, median),
                    integerResponse = FALSE
                )
                pdf(file = '../output/figures/DHARMa_northern_brms.beta.ry.disp.pdf')
                mod.resids %>% plot()
                dev.off()
                save(mod.resids, file=paste0('../data/modelled/resids.northern_brms.beta.ry.disp.RData'))
            }
            ## ----end
            dat.northern_brms.beta.ry.disp = emmeans(mod.northern_brms.beta.ry.disp,
                                                  ~Year, type='response') %>%
                as.data.frame()
            save(mod.northern_brms.beta.ry.disp, dat.northern_brms.beta.ry.disp, file=paste0('../data/modelled/mod.northern_brms.beta.ry.disp.RData'))
            rm(list=c('dat.northern_brms.beta.ry.disp', 'mod.northern_brms.beta.ry.disp'))
            gc()
        }  
    }
    ## ----end
    ## ---- Northern.MGCV.tow.beta
    {
        if ('MGCV beta' %in% models) {
            manta.tow.northern = manta.tow %>%
                filter(Region=='Northern GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL), nREEF_NAME=as.numeric(as.factor(REEF_NAME)))
            mod.northern_mgcv.beta <- gam(Cover ~ Year+s(nREEF_NAME, bs='re'),
                                          data=manta.tow.northern,
                                          family=betar,
                                          method='REML')

            dat.northern_mgcv.beta <- emmeans(mod.northern_mgcv.beta, ~Year, type='response') %>% as.data.frame() 
            save(dat.northern_mgcv.beta, mod.northern_mgcv.beta, file='../data/modelled/mod.northern_mgcv.beta.RData')
        }
    }
    ## ----end
    ## ---- Northern.glmmTMB.tow.beta vanilla
    {
        if ('glmmTMB_tow beta vanilla' %in% models) {
            mod.northern_glmmTMB.beta <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                                 data=manta.tow.northern,
                                                 ## weights=dat.all.northern$Tows,
                                                 family=beta_family())
            dat.northern_glmmTMB.beta = emmeans(mod.northern_glmmTMB.beta, ~Year, type='response') %>%
                as.data.frame()
            ## DHARMa::simulateResiduals(mod.northern_glmmTMB.beta, plot=TRUE)
            ## performance::check_model(mod.northern_glmmTMB.beta)
            save(mod.northern_glmmTMB.beta, dat.northern_glmmTMB.beta, file=paste0('../data/modelled/mod.northern_glmmTMB.beta.RData'))
        }
    }
    ## ----end
    ## ---- Northern.glmmTMB.tow.beta ry disp **
    {
        if ('glmmTMB_tow beta ry disp' %in% models & 'northern' %in% zone) {
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.northern_glmmTMB.beta.ry.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                                      dispformula = ~Year,
                                                      data=manta.tow.northern,
                                                      ## weights=dat.all.northern$Tows,
                                                      family=beta_family(),
                                                      control = glmmTMBControl(parallel = nt))
            dat.northern_glmmTMB.beta.ry.disp = emmeans(mod.northern_glmmTMB.beta.ry.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Northern.glmmTMB.tow.beta disp diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_northern_glmmTMB.beta.ry.disp.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.northern_glmmTMB.beta.ry.disp,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.northern_glmmTMB.beta.ry.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.northern_glmmTMB.beta.ry.disp.RData'))
            }
            ## ----end
            save(mod.northern_glmmTMB.beta.ry.disp, dat.northern_glmmTMB.beta.ry.disp,
                 file=paste0('../data/modelled/mod.northern_glmmTMB.beta.ry.disp.RData'))
        }
    }
    ## ----end
    ## ---- Northern.glmmTMB.tow.beta disp **
    {
        if ('glmmTMB_tow beta disp' %in% models & 'northern' %in% zone) {
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.northern_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME),
                                                      dispformula = ~Year,
                                                      data=manta.tow.northern,
                                                      ## weights=dat.all.northern$Tows,
                                                      family=beta_family(),
                                                      control = glmmTMBControl(parallel = nt))
            dat.northern_glmmTMB.beta.disp = emmeans(mod.northern_glmmTMB.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Northern.glmmTMB.tow.beta disp diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_northern_glmmTMB.beta.disp.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.northern_glmmTMB.beta.disp,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.northern_glmmTMB.beta.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.northern_glmmTMB.beta.disp.RData'))
            }
            ## ----end
            save(mod.northern_glmmTMB.beta.disp, dat.northern_glmmTMB.beta.disp,
                 file=paste0('../data/modelled/mod.northern_glmmTMB.beta.disp.RData'))
        }
    }
    ## ----end
    ## ---- Northern.glmmTMB.tow.beta disp random.effects
    {
        if ('glmmTMB_tow beta disp random slope' %in% models) {
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.northern_glmmTMB.beta.disp.rs <- glmmTMB(Cover ~ Year + (Year|REEF_NAME),
                                                         dispformula = ~Year,
                                                         data=manta.tow.northern,
                                                         ## weights=dat.all.northern$Tows,
                                                         family=beta_family(),
                                                         control = glmmTMBControl(parallel=nt))
            dat.northern_glmmTMB.beta.disp.rs = emmeans(mod.northern_glmmTMB.beta.disp.rs, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Northern.glmmTMB.tow.beta disp random.effects diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_northern_glmmTMB.beta.disp.rs.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.northern_glmmTMB.beta.disp.rs,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.northern_glmmTMB.beta.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.northern_glmmTMB.beta.disp.RData'))
                ## performance::check_model(mod.northern_glmmTMB.beta.disp.rs)
                ## ----end
            }
            save(mod.northern_glmmTMB.beta.disp.rs, dat.northern_glmmTMB.beta.disp.rs, file=paste0('../data/modelled/mod.northern_glmmTMB.beta.disp.rs.RData'))
        }
    }
    ## ----end
    ## ---- Northern.BRMS.tow.ordinal
    {
        if ('BRMS ordinal' %in% models) {
            library(brms)
            manta.tow.northern.tally <- manta.tow.northern %>%
                group_by(REEF_NAME, Year, oLIVE_CORAL) %>%
                count() %>%
                ungroup
            mod.northern_brms.cumulative <- brm(bf(oLIVE_CORAL|weights(n) ~ Year+(1|REEF_NAME)),
                                                data=manta.tow.northern.tally,
                                                ## family=cumulative("logit", threshold = 'flexible'),
                                                family=cumulative("logit", threshold = 'equidistant'),
                                                ## family=cumulative("probit", threshold = 'equidistant'),
                                                iter=1e4,
                                                warmup=5e3,
                                                thin=5,
                                                chains=4, cores=4,
                                                prior = prior(normal(0, 3), class = "b") +
                                                    prior(gamma(1, 0.5), class = "sd") +
                                                    prior(normal(0, 3), class = "Intercept") 
                                                )
            summary(mod.northern_brms.cumulative)
            ## save(mod.northern_brms.cumulative, file='../data/modelled/mod.northern_brms.cumulative.RData')

            ndata = manta.tow.northern %>% tidyr::expand(Year=levels(Year))
            pred1 <- posterior_epred(mod.northern_brms.cumulative, newdata=ndata, re_formula=NA) %>%
                aperm(c(1,3,2))
            lookup <-data.frame(LIVE_CORAL=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>%
                mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))
            ## lookup <- manta.tow.northern %>% tidyr::expand(LIVE_CORAL) %>%
            ##     mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))
            cats <- manta.tow.northern %>% group_by(LIVE_CORAL) %>% count() %>% ungroup() %>% pull(LIVE_CORAL)
            lookup1 = lookup %>% filter(LIVE_CORAL %in% cats) 
            out2=sweep(pred1, 2, lookup1$Cover, '*')
            out3 = apply(out2, 3, rowSums)
            dat.northern_brms.cumulative = ndata %>% cbind(tidyMCMC(as.mcmc(out3), conf.int=TRUE, conf.method='HPDinterval'))
            ## dat.northern_brms.cumulative %>%
            ##     ggplot(aes(y=estimate, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.3) +
            ##     geom_line() +
            ##     geom_ribbon(data=dat.northern, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.northern, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.northern.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.northern.original, aes(y=mean), color='red') +
            ##     geom_line(data=b, aes(y=Mean, x=as.numeric(as.character(Year))), color='purple', size=2)
            save(mod.northern_brms.cumulative, dat.northern_brms.cumulative, file='../data/modelled/mod.northern_brms.cumulative.RData')
        }  
    }
    ## ----end
    ## ---- Northern.MGCV.tow.ordinal
    {
        if ('MGCV ordinal' %in% models) {
            library(mgcv)
            manta.tow.northern = manta.tow %>%
                filter(Region=='Northern GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL), nREEF_NAME=as.numeric(as.factor(REEF_NAME)))
            mod.northern_mgcv.ordinal <- gam(nLIVE_CORAL ~ Year+s(nREEF_NAME, bs='re'),
                                             data=manta.tow.northern,
                                             family=ocat(R=16),
                                             method='REML')

            lookup <- manta.tow.northern %>% tidyr::expand(LIVE_CORAL) %>%
                mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))

            dat.northern_mgcv.ordinal <- predict(mod.northern_mgcv.ordinal, data.frame(Year=levels(manta.tow.northern$Year), nREEF_NAME=0), type='response', exclude='s(nREEF_NAME)') %>% as.data.frame %>%
                setNames(c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>% 
                mutate(Year=levels(manta.tow.northern$Year)) %>%
                pivot_longer(cols=-Year) %>%
                mutate(LIVE_CORAL=factor(name, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                         ordered=TRUE)) %>% 
                full_join(lookup) %>%
                mutate(P=value*Cover) %>%
                group_by(Year) %>%
                summarise(Mean=sum(P))

            ## dat.northern.mgcv %>%
            ##     ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            ##     ## geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), alpha=0.3) +
            ##     geom_line() +  
            ##     geom_ribbon(data=dat.northern, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.northern, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.northern.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.northern.original, aes(y=mean), color='red')
            save(mod.northern_mgcv.ordinal,dat.northern_mgcv.ordinal, file='../data/modelled/mod.northern_mgcv.ordinal.RData')
        }
    }
    ## ----end
    ## ---- Northern.CLMM.tow.ordinal
    {
        if ('CLMM' %in% models) {
            library(ordinal) 
            manta.tow.northern = manta.tow %>%
                filter(Region=='Northern GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE))
            mod.northern_clmm <-clmm(oLIVE_CORAL ~ Year + (1|REEF_NAME),
                                     data=manta.tow.northern)
            
            summary(mod.northern_clmm)
            dat.northern_clmm=emmeans(mod.northern_clmm, ~oLIVE_CORAL|Year, mode="prob") %>%
                as.data.frame() %>%
                group_by(Year) %>%
                mutate(oLIVE_CORAL=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>%
                ungroup %>%
                mutate(Perc=CoralTrends_calcPercent(oLIVE_CORAL),
                       Cover=prob*Perc) %>%
                group_by(Year) %>%
                summarise(Mean=sum(Cover))

            ## dat.northern_clmm %>% 
            ##     ggplot(aes(x=as.numeric(as.character(Year)))) +
            ##     geom_line(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(data=dat.northern, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.northern, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.northern.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.northern.original, aes(y=mean), color='red')

            save(mod.northern_clmm, dat.northern_clmm, file='../data/modelled/dat.northern_clmm.RData')
        }
    }
    ## ----end
    ## ---- Northern.glmmTMB.tow.beta.linear
    {
        if ('glmmTMB.beta.linear' %in% models) {
            mod.northern_glmmTMB.beta.linear <- glmmTMB(Cover ~ REPORT_YEAR + (1|REEF_NAME/REEF_YEAR),
                                                        data=manta.tow.northern,
                                                        ## weights=dat.all.northern$Tows,
                                                        family=beta_family())
            dat.northern_glmmTMB.beta.linear = emmeans(mod.northern_glmmTMB.beta.linear, ~REPORT_YEAR, at=list(REPORT_YEAR=unique(manta.tow.northern$REPORT_YEAR)), type='response') %>%
                as.data.frame()

            load(file=paste0('../data/modelled/mod.northern_glmmTMB.beta.disp.RData'))
            dat.northern_glmmTMB.beta.linear %>%
                ggplot() +
                geom_line(data=dat.northern_glmmTMB.beta.disp, aes(y=response, x=as.numeric(as.character(Year))), color='blue') +
                geom_ribbon(data=dat.northern_glmmTMB.beta.disp,aes(ymin=lower.CL, ymax=upper.CL, x=as.numeric(as.character(Year))), fill='lightblue', alpha=0.3) +
                geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, x=REPORT_YEAR), fill='orange', alpha=0.5) +
                geom_line(aes(y=response, x=REPORT_YEAR)) +
                scale_x_continuous('') +
                scale_y_continuous('Hard coral cover (%)', labels=function(x) x*100) +
                theme_classic() +
                ggtitle('Northern GBR')
        }
    }
    ## ----end

    ## Compare models
    {
        if (COMPARE_MODELS) {
            ## ---- beta.disp
        {
            load(file=paste0('../data/modelled/mod.northern_glmmTMB.beta.disp.RData'))
            load(file=paste0('../data/modelled/mod.northern_glmmTMB.beta.disp.rs.RData'))
            load(file=paste0('../data/modelled/mod.northern_brms.beta.disp.RData'))
            load(file=paste0('../data/modelled/mod.northern_inla.beta.disp.RData'))

            pdf(file = '../output/figures/comparison.northern.pdf')
            g1 <- ggplot() +
                geom_line(data = dat.northern_brms.beta.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'brms')) +
                geom_ribbon(data = dat.northern_brms.beta.disp,
                            aes(y = response,
                                ymin = lower.HPD, ymax = upper.HPD,
                                x = as.numeric(as.character(Year)),
                                fill = 'brms'),
                            alpha=0.3) +
                geom_line(data = dat.northern_glmmTMB.beta.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'glmmTMB')) +
                geom_ribbon(data = dat.northern_glmmTMB.beta.disp,
                            aes(y = response,
                                ymin = lower.CL, ymax = upper.CL,
                                x = as.numeric(as.character(Year)),
                                fill = 'glmmTMB'),
                            alpha=0.3) + 
                geom_line(data = dat.northern_inla.beta.disp,
                          aes(y = mean,
                              x = as.numeric(as.character(Year)),
                              colour = 'inla')) +
                geom_ribbon(data = dat.northern_inla.beta.disp,
                            aes(y = mean,
                                ymin = lower, ymax = upper,
                                x = as.numeric(as.character(Year)),
                                fill = 'inla'),
                            alpha=0.3) 
            print(g1)
            dev.off()
        }
            ## ----end
            ## ---- beta.ry.disp
        {
            load(file=paste0('../data/modelled/mod.northern_glmmTMB.beta.ry.disp.RData'))
            load(file=paste0('../data/modelled/mod.northern_brms.beta.ry.disp.RData'))
            load(file=paste0('../data/modelled/mod.northern_inla.beta.ry.disp.RData'))
            pdf(file = '../output/figures/comparison.northern.ry.pdf')
            g1 <- ggplot() +
                geom_line(data = dat.northern_brms.beta.ry.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'brms')) +
                geom_ribbon(data = dat.northern_brms.beta.ry.disp,
                            aes(y = response,
                                ymin = lower.HPD, ymax = upper.HPD,
                                x = as.numeric(as.character(Year)),
                                fill = 'brms'),
                            alpha=0.3) +
                geom_line(data = dat.northern_glmmTMB.beta.ry.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'glmmTMB')) +
                geom_ribbon(data = dat.northern_glmmTMB.beta.ry.disp,
                            aes(y = response,
                                ymin = lower.CL, ymax = upper.CL,
                                x = as.numeric(as.character(Year)),
                                fill = 'glmmTMB'),
                            alpha=0.3) + 
                geom_line(data = dat.northern_inla.beta.ry.disp,
                          aes(y = mean,
                              x = as.numeric(as.character(Year)),
                              colour = 'inla')) +
                geom_ribbon(data = dat.northern_inla.beta.ry.disp,
                            aes(y = mean,
                                ymin = lower, ymax = upper,
                                x = as.numeric(as.character(Year)),
                                fill = 'inla'),
                            alpha=0.3) 
            print(g1)
            dev.off()
        }
            ## ----end

        }
    }
    
    if (1==2) {
        ## Compare the models
    {
        ## ---- Northern Compare models

        load(file='../data/modelled/dat.northern.RData')
        dat.northern.original <- dat.northern
        load(file='../data/modelled/mod.northern.RData')
        load(file='../data/modelled/mod.northern_inla_beta.RData')
        load(file='../data/modelled/newdata.northern_inla_beta.RData')
        newdata.northern_beta <- newdata.northern
        load(file='../data/modelled/mod.northern_inla_binomial.RData')
        load(file='../data/modelled/newdata.northern_inla_binomial.RData')
        newdata.northern_binomial <- newdata.northern
        load(file='../data/modelled/dat.northern_inla_tow.RData')
        load(file='../data/modelled/dat.northern_mgcv.RData')
        load(file='../data/modelled/dat.northern_brms.cumulative.RData')
        load(file='../data/modelled/dat.northern_clmm.RData')
        load(file='../data/modelled/dat.northern_brms.beta.RData')
        load(file='../data/modelled/dat.northern_glmmTMB.RData')
        load(file='../data/modelled/dat.northern_mgcv.beta.RData')

        ## original
        g1 <- dat.northern.original %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('Original (stan binomial reef level)') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g1

        ## inla (reef level binomial)
        g2 <- newdata.northern_binomial %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA reef level binomial') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g2

        ## inla (reef level beta)
        g3 <- newdata.northern_beta %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA reef level beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g3

        ## inla (tow level beta)
        g4 <- dat.northern %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA tow level beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g4

        ## mgcv ordinal
        g5 <- dat.northern.mgcv %>%
            ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('mgcv ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))

        ## brms ordinal
        g6 <- dat.northern_brms.cumulative %>%
            ggplot(aes(y=estimate, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill='blue', alpha=0.3) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g6    

        ## clmm ordinal
        g7 <- dat.northern_clmm %>%
            ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('clmm ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g7    

        ## brms beta
        g8 <- dat.northern_brms.beta %>%
            ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g8

        library(patchwork)
        g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8


        dat.northern_brms.beta %>%
            ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ## geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill='BRMS beta'), alpha=0.3) +
            geom_line(aes(color='BRMS beta')) +
            ## geom_ribbon(data=dat.northern_brms.cumulative, aes(y=estimate, ymin=conf.low, ymax=conf.high, fill='BRMS ordinal'), alpha=0.3) +
            geom_line(data=dat.northern_brms.cumulative, aes(y=estimate, color='BRMS ordinal')) +
            ## geom_ribbon(data=newdata.northern_beta, aes(y=mean, ymin=lower, ymax=upper, fill='INLA beta'), alpha=0.3) +
            geom_line(data=newdata.northern_beta, aes(y=mean, color='INLA beta')) +
            geom_line(data=dat.northern_clmm, aes(y=Mean, color='clmm ordinal')) +
            geom_line(data=dat.northern.mgcv, aes(y=Mean, color='mgcv ordinal')) +
            geom_line(data=dat.northern_glmmTMB, aes(y=response, color='glmmTMB beta')) +
            geom_line(data=dat.northern_mgcv.beta, aes(y=response, color='mgcv beta')) +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_brewer('Model', type='qual') + 
            scale_fill_discrete('Model') + 
            ## rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))


        rm(list=c('dat.northern','mod.northern','mod.northern_inla_beta','newdata.northern','newdata.northern_beta', 'mod_northern_inla_binomial','newdata.northern_binomial'))
    }
        ## ----end
        ## ---- junk
    {

        ## summary(mod.northern)
        ## dat.northern = data.frame(Location='Northern',Year=unique(dat.all.northern$Year), N=length(unique(dat.all.northern$REEF_NAME)))
        ## Xmat = model.matrix(~Year, dat.northern)
        ## coefs = data.frame(mod.northern) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
        ## coefs = mod.northern$fit %>% as.data.frame() %>% dplyr::select(matches('b_.*'))
        ## Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
        ## dat.northern = cbind(dat.northern,
        ##             plyr:::adply(Fit,2,function(x) {
        ##     data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
        ##             })
        ## )


        ##     ## glmmTMB
        ##     library(glmmTMB)
        ##     library(emmeans)
        ##     mod.glmmTMB <-  glmmTMB(Cover ~ Year+(Year|REEF_NAME),
        ##                        data=dat.all.northern,
        ##                        #weights=dat.all.gbr$Weight,
        ##                        family=beta_family)
        ##     emmeans(mod.glmmTMB, ~Year, type='response') %>% as.data.frame %>%
        ##     ggplot() +
        ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, x=as.numeric(as.character(Year))), fill='blue', alpha=0.2) +
        ##     geom_line(aes(y=response, x=as.numeric(as.character(Year)))) +
        ##     theme_bw()


        ## dat.northern %>% arrange(Year)
        ## dat.northern %>%
        ##     ggplot(aes(y=mean, x=as.numeric(Year))) +
        ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        ##     geom_line()

        ## save(dat.northern, file='../data/modelled/dat.northern.RData')
        ## save(dat.all.northern, file='../data/modelled/dat.all.northern.RData')
        ## save(mod.northern, file='../data/modelled/mod.northern.RData')
        ## rm(list=c('last_year','mod.northern', 'dat.northern', 'coefs', 'Fit'))
        ## gc()
    }
        ## ----end
    }
    ## ----end
}

    
## Central
{
    ## ---- Central
    ## ---- Central.Data
    {
        dat.all.central = dat.all %>%
            filter(Location=='Central GBR') %>%
            droplevels %>% 
            mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
            group_by(REEF_NAME) %>%
            mutate(W=mean(Tows, na.rm=TRUE)) %>%
            ungroup %>%
            mutate(W1=W/sum(W)) %>%
            group_by(Year) %>%
            mutate(W2=Tows/sum(Tows)) %>%
            ungroup
        save(dat.all.central, file='../data/modelled/dat.all.central.RData')

        ## Tow level data
        manta.tow.central = manta.tow %>%
            filter(Region=='Central GBR') %>%
            droplevels %>%
            mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                                      levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                      ordered=TRUE),
                   nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
                   nLIVE_CORAL=as.numeric(oLIVE_CORAL),
                   REEF_YEAR = interaction(REEF_NAME, Year)
                   )
        save(manta.tow.central, file='../data/modelled/manta.tow.central.RData')
    }
    ## ----end
    ## ---- Central.Raw cells
    {
        dat.all.central.cellmeans <- cellMeansRaw(dat.all.central)
        rawAdd <- ggproto_Raw(dat.all.central.cellmeans)

        manta.tow.central.cellmeans <- cellMeansRaw(manta.tow %>%
                                                    filter(Region=='Central GBR') %>%
                                                    droplevels %>% 
                                                    group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                                    summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                                    ungroup)
    }
    ## ----end                                        

    ## ---- Central.stan_glmer.reef.binomial
    {
        if ('original' %in% models) {
            ## mod.central <- ModelOriginal(form=formula(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME)),
            ##                               dat=dat.all.central,
            ##                               location='Central')
            ## mod.central[['newdata']] %>%
            ##     ggplot(aes(y=mean, x=as.numeric(Year))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line()
            ## mod.central[['mod']] %>% save(file='../data/modelled/mod.central_original.RData')
            ## mod.central[['newdata']] %>% save(file='../data/modelled/newdata.central_original.RData')

            mod.central_glmer <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                             data=dat.all.central, family=binomial,iter=5000,warmup=2500,
                                             chains=3,cores=3)
            dat.central_glmer <- data.frame(Location='Central GBR',Year=unique(dat.all.central$Year), N=length(unique(dat.all.central$REEF_NAME)))
            Xmat <- model.matrix(~Year, dat.central_glmer)

            coefs = data.frame(mod.central_glmer) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
            Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
            dat.central_glmer = cbind(dat.central_glmer,
                                      plyr:::adply(Fit,2,function(x) {
                                          data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                                      })
                                      )
            save(dat.central_glmer, mod.central_glmer, file='../data/modelled/dat.central_glmer.RData')
            ## save(dat.central, file='../data/modelled/dat.central.RData')
            ## save(mod.central, file='../data/modelled/mod.central.RData')
            rm(list=c('last_year','dat.central_glmer','mod.central_glmer'))
            gc()
            ##   mod.central <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
            ##                          data=dat.all.central, family=binomial,iter=5000,warmup=2500,
            ##                          chains=3,cores=3)
            ##   dat.central <- data.frame(Location='Central GBR',Year=unique(dat.all.central$Year), N=length(unique(dat.all.central$REEF_NAME)))
            ##   Xmat <- model.matrix(~Year, dat.central)

            ##   coefs = data.frame(mod.central) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
            ##   Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
            ##   dat.central = cbind(dat.central,
            ##                   plyr:::adply(Fit,2,function(x) {
            ##                       data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
            ##                   })
            ##                   )
            ##   save(dat.central, file='../data/modelled/dat.central.RData')
            ##   save(mod.central, file='../data/modelled/mod.central.RData')
            ##   rm(list=c('last_year','mod.central'))
            ##   gc()
            ## ## mod.central <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
            ## ##                             data=dat.all.central,
            ## ##                             family=binomial,
            ## ##                             iter=5000,
            ## ##                             warmup=2500,
            ##   ##                             chains=3,cores=3)
        }
    }
    ## ----end
    ## ---- Central.stan_glmer.reef.beta
    {
        if ('stan_glmer beta' %in% models) {
            mod.central <- Modelstan_glmer(form=Cover ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                           dat=dat.all.central)
            newdata.central <- cellMeansOriginal(mod.central, dat=dat.all.central, location='Central') 
            newdata.central %>%
                ggplot(aes(y=mean, x=as.numeric(Year))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                geom_line()
            save(mod.central, file='../data/modelled/mod.central_stan_glmer.RData')
            save(newdata.central, file='../data/modelled/newdata.central_stan_glmer.RData')
            rm(list=c('last_year','mod.central', 'newdata.central', 'coefs', 'Fit'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.BRMS.reef.beta
    {
        if ('BRMS_reef beta' %in% models) {
            priors <- prior(normal(0, 5), class = "b") +
                prior(normal(0, 5), class = "Intercept") +
                ## prior(gamma(2, 0.1), class = "sd") +
                prior(gamma(1, 0.5), class = "sd") +
                prior(gamma(0.01, 0.01), class = "phi")
            inits = list(list(phi=list(rgamma(1,0.1,0.1))),
                         list(phi=list(rgamma(1,0.1,0.1))),
                         list(phi=list(rgamma(1,0.1,0.1)))
                         )
            mod.central <- Modelbrms(form=Cover|weights(Tows) ~ Year + (1|REEF_NAME),
                                     dat=dat.all.central)
            
            mod.central <- Modelbrms(form=Cover|weights(Tows) ~ Year + (Year|REEF_NAME),
                                     dat=dat.all.central)
        }
    }
    ## ----end
    ## ---- Central.INLA.reef.beta vanilla
    {
        if ("INLA_reef beta" %in% models) {
            dat.inla <- dataINLA(dat=dat.all.central %>% mutate(W=NA))
            dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
            mod.central_inla_reef.beta <- ModelINLA_beta(form=Cover~Year+
                                                             f(P_CODE.mod, model='iid')+
                                                             f(REEF_NAME, model='iid'),
                                                         dat=dat.inla[['dat.1']])
            dat.central_inla_reef.beta <- cellMeansINLA(mod=mod.central_inla_reef.beta, newdata.hcc=dat.inla[['newdata.hcc']],
                                                        n.2=dat.inla[['n.2']])
            ## newdata.central %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd
            save(dat.central_inla_reef.beta, mod.central_inla_reef.beta, file='../data/modelled/mod.central_inla_reef.beta.RData')
            ## save(newdata.central, file='../data/modelled/newdata.central_inla_beta.RData')
            rm(list=c('mod.central_inla_reef.beta', 'dat.central_inla_reef.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.INLA.reef.beta scaled
    {
        if ("INLA_reef beta" %in% models) {
            dat.inla <- dataINLA(dat=dat.all.central %>% mutate(W=NA))
            dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
            mod.central_inla_reef.beta.scaled <- ModelINLA_beta(form=Cover~Year+
                                                                    f(P_CODE.mod, model='iid')+
                                                                    f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                                                    f(REEF_NAME1, Year, model='iid', scale=dat.scale$Tows),
                                                                dat=dat.inla[['dat.1']])
            dat.central_inla_reef.beta.scaled <- cellMeansINLA(mod=mod.central_inla_reef.beta.scaled, newdata.hcc=dat.inla[['newdata.hcc']],
                                                               n.2=dat.inla[['n.2']])
            ## newdata.central %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd
            save(dat.central_inla_reef.beta.scaled, mod.central_inla_reef.beta.scaled, file='../data/modelled/mod.central_inla_reef.beta.scaled.RData')
            ## save(newdata.central, file='../data/modelled/newdata.central_inla_beta.RData')
            rm(list=c('mod.central_inla_reef.beta.scaled', 'dat.central_inla_reef.beta.scaled'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.INLA.reef.binomial
    {
        if ('INLA_reef binomial' %in% models) {
            dat.inla <- dataINLA(dat=dat.all.central)
            mod.central <- ModelINLA_binomial(form=Cvr1~Year+
                                                  f(P_CODE.mod, model='iid')+
                                                  f(REEF_NAME, model='iid') +
                                                  f(REEF_NAME1, Year, model='iid'),
                                              dat=dat.inla[['dat.1']])
            newdata.central <- cellMeansINLA(mod=mod.central, newdata.hcc=dat.inla[['newdata.hcc']],
                                             n.2=dat.inla[['n.2']])
            newdata.central %>%
                ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                geom_line() +
                rawAdd
            save(mod.central, file='../data/modelled/mod.central_inla_binomial.RData')
            save(newdata.central,file='../data/modelled/newdata.central_inla_binomial.RData')
            rm(list=c('mod.central', 'newdata.central'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.glmmTMB.reef.beta
    {
        if ('glmmTMB_reef beta' %in% models) {
            mod.central_glmmTMB.reef.beta <- glmmTMB(Cover ~ Year + (1|P_CODE.mod) +  (1|REEF_NAME),
                                                     data=dat.all.central,
                                                     weights=dat.all.central$Tows,
                                                     family=beta_family())
            dat.central_glmmTMB.reef.beta  = emmeans(mod.central_glmmTMB.reef.beta, ~Year, type='response') %>%
                as.data.frame()
            save(mod.central_glmmTMB.reef.beta, dat.central_glmmTMB.reef.beta, file='../data/modelled/dat.central_glmmTMB.reef.beta.RData')
            
            ## emmeans(mod.central, ~Year, type='response') %>%
            ##     as.data.frame() %>% 
            ##     ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd

            rm(list=c('dat.central_glmmTMB.reef.beta', 'mod.central_glmmTMB.reef.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.glmmTMB.reef.beta disp
    {
        if ('glmmTMB_reef beta' %in% models) {
            mod.central_glmmTMB.reef.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME),
                                                          dispformula = ~Year,
                                                          data=dat.all.central,
                                                          weights=dat.all.central$Tows,
                                                          family=beta_family())
            dat.central_glmmTMB.reef.beta.disp  = emmeans(mod.central_glmmTMB.reef.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            save(mod.central_glmmTMB.reef.beta.disp, dat.central_glmmTMB.reef.beta.disp, file='../data/modelled/dat.central_glmmTMB.reef.beta.disp.RData')
            
            ## emmeans(mod.central, ~Year, type='response') %>%
            ##     as.data.frame() %>% 
            ##     ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd

            rm(list=c('dat.central_glmmTMB.reef.beta.disp', 'mod.central_glmmTMB.reef.beta.disp'))
            gc()
        }
    }
    ## ----end

    ## ---- Central.INLA.tow.beta scaled
    {
        if ("INLA_tow beta scaled" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.central, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.central %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            mod.central_inla.beta.scaled <- ModelINLA_beta(form=Cover~Year +
                                                               f(P_CODE.mod, model='iid') +
                                                               f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                                               f(REEF_NAME2, YEAR1, model='iid', scale=dat.scale$Tows),
                                                           dat=dat,
                                                           family='beta',
                                                           weights=NULL)
            
            dat.central_inla.beta.scaled <- cellMeansINLA(mod=mod.central_inla.beta.scaled, newdata.hcc=dat.inla[['newdata.hcc']],
                                                          n.2=dat.inla[['n.2']], FUN=plogis)
            ## dat.raw = dat %>% group_by(REEF_NAME, Year) %>% summarize(Cover=mean(Cover, na.rm=TRUE)) %>% filter(!is.na(Cover)) %>% droplevels
            ## dat.central %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill='blue') +
            ##     geom_line(color='blue') +
            ##    rawAdd +
            ##    ggtitle('GBR INLA tow level') #
            ##    ## geom_line(data=mod.brms.df, aes(y=response, x=as.numeric(as.character(year)))) +
            ##    ## geom_ribbon(data=mod.brms.df, aes(y=response, ymin=lower.HPD, ymax=upper.HPD, x=as.numeric(as.character(year))), alpha=0.3)        

            save(mod.central_inla.beta.scaled, dat.central_inla.beta.scaled, file='../data/modelled/mod.central_inla.beta.scaled.RData')
            rm(list=c('dat.central_inla.beta.scaled', 'mod.central_inla.beta.scaled'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.INLA.tow.beta
    {
        if ("INLA_tow beta" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.central, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.central %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            mod.central_inla.beta <- ModelINLA_beta(form=Cover~Year +
                                                        f(P_CODE.mod, model='iid') +
                                                        f(REEF_NAME, model='iid') + 
                                                        f(REEF_YEAR, model='iid'), 
                                                    dat=dat,
                                                    family='beta',
                                                    weights=NULL)
            
            dat.central_inla.beta <- cellMeansINLA(mod=mod.central_inla.beta, newdata.hcc=dat.inla[['newdata.hcc']],
                                                   n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.central_inla.beta, dat.central_inla.beta, file='../data/modelled/mod.central_inla.beta.RData')
            rm(list=c('dat.central_inla.beta', 'mod.central_inla.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.INLA.tow.beta **
    {
        if ("INLA_tow beta disp" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.central, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.central %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            ## we will leverage a mixed likelihood model
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            
            dd <- dat %>%
                dplyr::select(Cover, Year, REEF_NAME) %>%
                mutate(YEAR = Year)%>%
                pivot_wider(id_cols = c(YEAR,REEF_NAME),
                            names_from = Year,
                            values_from = Cover)
            dd1 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(matches("[0-9]{4}")) %>%
                as.matrix()
            dd2 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(YEAR, REEF_NAME)

            mod.central_inla.beta.disp <- inla(form = dd1~YEAR +
                                               f(REEF_NAME, model='iid'),
                                           dat=dd2,
                                           family=rep('beta',ncol(dd1)),
                                           control.fixed = list(mean = 0, prec = 0.001,
                                                                mean.intercept = 0.5,
                                                                prec.intercept = 0.001),
                                           control.predictor = list(compute = TRUE,
                                                                    link = 1,
                                                                    quantiles = c(0.025,0.25,0.5,0.75,0.975)
                                                                    )
                                           )
            dat.central_inla.beta.disp <- cellMeansINLA(mod=mod.central_inla.beta.disp, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.central_inla.beta.disp, dat.central_inla.beta.disp, file='../data/modelled/mod.central_inla.beta.disp.RData')
            rm(list=c('dat.central_inla.beta.disp', 'mod.central_inla.beta.disp'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.INLA.tow.ry.beta **
    {
        if ("INLA_tow beta ry disp" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.central, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.central %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            ## we will leverage a mixed likelihood model
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            
            dd <- dat %>%
                dplyr::select(Cover, Year, REEF_NAME, REEF_YEAR) %>%
                mutate(YEAR = Year)%>%
                pivot_wider(id_cols = c(YEAR,REEF_NAME, REEF_YEAR),
                            names_from = Year,
                            values_from = Cover)
            dd1 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(matches("[0-9]{4}")) %>%
                as.matrix()
            dd2 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(YEAR, REEF_NAME, REEF_YEAR)

            mod.central_inla.beta.ry.disp <- inla(form = dd1~YEAR +
                                               f(REEF_NAME, model='iid') +
                                               f(REEF_YEAR, model='iid'),
                                           dat=dd2,
                                           family=rep('beta',ncol(dd1)),
                                           control.fixed = list(mean = 0, prec = 0.001,
                                                                mean.intercept = 0.5,
                                                                prec.intercept = 0.001),
                                           control.predictor = list(compute = TRUE,
                                                                    link = 1,
                                                                    quantiles = c(0.025,0.25,0.5,0.75,0.975)
                                                                    )
                                           )
            dat.central_inla.beta.ry.disp <- cellMeansINLA(mod=mod.central_inla.beta.ry.disp, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.central_inla.beta.ry.disp, dat.central_inla.beta.ry.disp, file='../data/modelled/mod.central_inla.beta.ry.disp.RData')
            rm(list=c('dat.central_inla.beta.ry.disp', 'mod.central_inla.beta.ry.disp'))
            gc()
        }
    }
    ## ----end
    ## ---- Central.BRMS.tow.beta vanilla
    {
        if ('BRMS beta vanilla' %in% models) {
            mod.central_brms.beta <- brm(Cover ~ Year + (1|REEF_NAME),
                                         data=manta.tow.central,
                                         family=Beta(link='logit'),
                                         iter=1e4,
                                         warmup=5e3,
                                         thin=5,
                                         chains=4, cores=4,
                                         prior = prior(normal(0, 3), class = "b") +
                                             prior(normal(0, 3), class = "Intercept") +
                                             prior(gamma(2, 1), class = "sd") +
                                             prior(gamma(2, 1), class = "sd", group = "REEF_NAME") +
                                             prior(gamma(2, 1), class = "phi")
                                         )
            dat.central_brms.beta = emmeans(mod.central_brms.beta, ~Year, type='response') %>%
                as.data.frame()
            save(mod.central_brms.beta, dat.central_brms.beta, file=paste0('../data/modelled/mod.central_brms_beta.RData'))
            ## dat.central_brms.beta %>%
            ##     ggplot(aes(y=response, x=as.numeric(as.character(year)))) +
            ##     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), alpha=0.3) +
            ##     geom_line() +
            ##     geom_ribbon(data=dat.central, aes(y=mean, x=as.numeric(as.character(Year)), ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.central, aes(y=mean, x=as.numeric(as.character(Year))), color='blue') +
            ##     geom_ribbon(data=dat.central.original, aes(y=mean, x=as.numeric(as.character(Year)), ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.central.original, aes(y=mean, x=as.numeric(as.character(Year))), color='red')
            ## save(dat.central_brms.beta, file='../data/modelled/dat.central_brms.beta.RData')
            rm(list=c('dat.central_brms.beta', 'mod.central_brms.beta'))
            gc()
        }    
    }
    ## ----end
    ## ---- Central.BRMS.tow.beta disp **
    {
        if ('BRMS beta disp' %in% models & 'central' %in% zone) {
            cat('Fitting brms disp for Central\n\n')
            mod.central_brms.beta.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME), phi~0+Year),
                                              data=manta.tow.central,
                                              family=Beta(link='logit'),
                                              iter=1e4,
                                              warmup=5e3,
                                              thin=5,
                                              chains=4, cores=4,
                                              prior = prior(normal(0, 3), class = "b") +
                                                  prior(normal(0, 3), class = "Intercept") +
                                                  prior(gamma(2, 1), class = "sd") #+
                                              ## prior(gamma(2, 1), class = "phi")
                                              )
            dat.central_brms.beta.disp = emmeans(mod.central_brms.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Central.BRMS.tow.beta disp diagnostics
            {
                ## sampling diagnostics
                pdf(file = '../output/figures/traceplots_central_brms.beta.disp.pdf')
                rstan::traceplot(mod.central_brms.beta.disp$fit)
                dev.off()

                ## density overlay
                pdf(file = '../output/figures/density_central_brms.beta.disp.pdf')
                mod.central_brms.beta.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100) 
                dev.off()

                ## DHARMa residuals
                preds <- mod.central_brms.beta.disp %>%
                    posterior_predict(ndraws = 250, summary = FALSE)
                mod.resids <- createDHARMa(
                    simulatedResponse = t(preds),
                    observedResponse = manta.tow.central$Cover,
                    fittedPredictedResponse = apply(preds, 2, median),
                    integerResponse = FALSE
                )
                pdf(file = '../output/figures/DHARMa_central_brms.beta.disp.pdf')
                mod.resids %>% plot()
                dev.off()
                save(mod.resids, file=paste0('../data/modelled/resids.central_brms.beta.disp.RData'))
            }
            ## ----end
            save(mod.central_brms.beta.disp, dat.central_brms.beta.disp, file=paste0('../data/modelled/mod.central_brms.beta.disp.RData'))
            rm(list=c('dat.central_brms.beta.disp', 'mod.central_brms.beta.disp'))
            gc()
        }  
    }
    ## ----end
    ## ---- Central.BRMS.tow.beta ry disp **
    {
        if ('BRMS beta ry disp' %in% models & 'central' %in% zone) {
            cat('Fitting brms ry disp for Central\n\n')
            mod.central_brms.beta.ry.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME/REEF_YEAR), phi~0+Year),
                                              data=manta.tow.central,
                                              family=Beta(link='logit'),
                                              iter=1e4,
                                              warmup=5e3,
                                              thin=5,
                                              chains=4, cores=4,
                                              prior = prior(normal(0, 3), class = "b") +
                                                  prior(normal(0, 3), class = "Intercept") +
                                                  prior(gamma(2, 1), class = "sd") #+
                                              ## prior(gamma(2, 1), class = "phi")
                                              )
            dat.central_brms.beta.ry.disp = emmeans(mod.central_brms.beta.ry.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Central.BRMS.tow.beta ry disp diagnostics
            {
                ## sampling diagnostics
                pdf(file = '../output/figures/traceplots_central_brms.beta.ry.disp.pdf')
                rstan::traceplot(mod.central_brms.beta.ry.disp$fit)
                dev.off()

                ## density overlay
                pdf(file = '../output/figures/density_central_brms.beta.ry.disp.pdf')
                mod.central_brms.beta.ry.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100) 
                dev.off()

                ## DHARMa residuals
                preds <- mod.central_brms.beta.ry.disp %>%
                    posterior_predict(ndraws = 250, summary = FALSE)
                mod.resids <- createDHARMa(
                    simulatedResponse = t(preds),
                    observedResponse = manta.tow.central$Cover,
                    fittedPredictedResponse = apply(preds, 2, median),
                    integerResponse = FALSE
                )
                pdf(file = '../output/figures/DHARMa_central_brms.beta.ry.disp.pdf')
                mod.resids %>% plot()
                dev.off()
                save(mod.resids, file=paste0('../data/modelled/resids.central_brms.beta.ry.disp.RData'))
            }
            ## ----end
            save(mod.central_brms.beta.ry.disp, dat.central_brms.beta.ry.disp, file=paste0('../data/modelled/mod.central_brms.beta.ry.disp.RData'))
            rm(list=c('dat.central_brms.beta.ry.disp', 'mod.central_brms.beta.ry.disp'))
            gc()
        }  
    }
    ## ----end
    ## ---- Central.MGCV.tow.beta
    {
        if ('MGCV beta' %in% models) {
            manta.tow.central = manta.tow %>%
                filter(Region=='Central GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL), nREEF_NAME=as.numeric(as.factor(REEF_NAME)))
            mod.central_mgcv.beta <- gam(Cover ~ Year+s(nREEF_NAME, bs='re'),
                                         data=manta.tow.central,
                                         family=betar,
                                         method='REML')

            dat.central_mgcv.beta <- emmeans(mod.central_mgcv.beta, ~Year, type='response') %>% as.data.frame() 
            save(dat.central_mgcv.beta, mod.central_mgcv.beta, file='../data/modelled/mod.central_mgcv.beta.RData')
        }
    }
    ## ----end
    ## ---- Central.glmmTMB.tow.beta vanilla
    {
        if ('glmmTMB_tow beta vanilla' %in% models) {
            mod.central_glmmTMB.beta <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                                data=manta.tow.central,
                                                ## weights=dat.all.central$Tows,
                                                family=beta_family())
            dat.central_glmmTMB.beta = emmeans(mod.central_glmmTMB.beta, ~Year, type='response') %>%
                as.data.frame()
            ## DHARMa::simulateResiduals(mod.central_glmmTMB.beta, plot=TRUE)
            ## performance::check_model(mod.central_glmmTMB.beta)
            save(mod.central_glmmTMB.beta, dat.central_glmmTMB.beta, file=paste0('../data/modelled/mod.central_glmmTMB.beta.RData'))
        }
    }
    ## ----end
    ## ---- Central.glmmTMB.tow.beta disp **
    {
        if ('glmmTMB_tow beta disp' %in% models & 'central' %in% zone) {
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.central_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME),
                                                     dispformula = ~Year,
                                                     data=manta.tow.central,
                                                     ## weights=dat.all.central$Tows,
                                                     family=beta_family(),
                                                     control = glmmTMBControl(parallel = nt))
            dat.central_glmmTMB.beta.disp = emmeans(mod.central_glmmTMB.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Central.glmmTMB.tow.beta disp diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_central_glmmTMB.beta.disp.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.central_glmmTMB.beta.disp,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.central_glmmTMB.beta.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.central_glmmTMB.beta.disp.RData'))
            }
            ## ----end
            save(mod.central_glmmTMB.beta.disp, dat.central_glmmTMB.beta.disp, file=paste0('../data/modelled/mod.central_glmmTMB.beta.disp.RData'))
        }
    }
    ## ----end
    ## ---- Central.glmmTMB.tow.beta ry disp **
    {
        if ('glmmTMB_tow beta ry disp' %in% models & 'central' %in% zone) {
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.central_glmmTMB.beta.ry.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                                     dispformula = ~Year,
                                                     data=manta.tow.central,
                                                     ## weights=dat.all.central$Tows,
                                                     family=beta_family(),
                                                     control = glmmTMBControl(parallel = nt))
            dat.central_glmmTMB.beta.ry.disp = emmeans(mod.central_glmmTMB.beta.ry.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Central.glmmTMB.tow.beta ry disp diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_central_glmmTMB.beta.ry.disp.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.central_glmmTMB.beta.ry.disp,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.central_glmmTMB.beta.ry.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.central_glmmTMB.beta.ry.disp.RData'))
            }
            ## ----end
            save(mod.central_glmmTMB.beta.ry.disp, dat.central_glmmTMB.beta.ry.disp, file=paste0('../data/modelled/mod.central_glmmTMB.beta.ry.disp.RData'))
        }
    }
    ## ----end
    ## ---- Central.glmmTMB.tow.beta disp random.effects
    {
        if ('glmmTMB_tow beta disp' %in% models) {
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.central_glmmTMB.beta.disp.rs <- glmmTMB(Cover ~ Year + (Year|REEF_NAME),
                                                        dispformula = ~Year,
                                                        data=manta.tow.central,
                                                        ## weights=dat.all.central$Tows,
                                                        family=beta_family(),
                                                        control = glmmTMBControl(parallel=nt))
            dat.central_glmmTMB.beta.disp.rs = emmeans(mod.central_glmmTMB.beta.disp.rs,
                                                       ~Year, type='response') %>%
                as.data.frame()
            ## did not converge
            ## ---- Central.glmmTMB.tow.beta disp random.effects diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_central_glmmTMB.beta.disp.rs.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.central_glmmTMB.beta.disp.rs,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.central_glmmTMB.beta.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.central_glmmTMB.beta.disp.RData'))
                ## performance::check_model(mod.central_glmmTMB.beta.disp.rs)
            }
            ## ----end
            save(mod.central_glmmTMB.beta.disp.rs, dat.central_glmmTMB.beta.disp.rs, file=paste0('../data/modelled/mod.central_glmmTMB.beta.disp.rs.RData'))
        }
    }
    ## ----end
    ## ---- Central.BRMS.tow.ordinal
    {
        if ('BRMS ordinal' %in% models) {
            library(brms)
            manta.tow.central.tally <- manta.tow.central %>%
                group_by(REEF_NAME, Year, oLIVE_CORAL) %>%
                count() %>%
                ungroup
            mod.central_brms.cumulative <- brm(bf(oLIVE_CORAL|weights(n) ~ Year+(1|REEF_NAME)),
                                               data=manta.tow.central.tally,
                                               ## family=cumulative("logit", threshold = 'flexible'),
                                               family=cumulative("logit", threshold = 'equidistant'),
                                               ## family=cumulative("probit", threshold = 'equidistant'),
                                               iter=1e4,
                                               warmup=5e3,
                                               thin=5,
                                               chains=4, cores=4,
                                               prior = prior(normal(0, 3), class = "b") +
                                                   prior(gamma(1, 0.5), class = "sd") +
                                                   prior(normal(0, 3), class = "Intercept") 
                                               )
            summary(mod.central_brms.cumulative)
            ## save(mod.central_brms.cumulative, file='../data/modelled/mod.central_brms.cumulative.RData')

            ndata = manta.tow.central %>% tidyr::expand(Year=levels(Year))
            pred1 <- posterior_epred(mod.central_brms.cumulative, newdata=ndata, re_formula=NA) %>%
                aperm(c(1,3,2))
            lookup <-data.frame(LIVE_CORAL=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>%
                mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))
            ## lookup <- manta.tow.central %>% tidyr::expand(LIVE_CORAL) %>%
            ##     mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))
            cats <- manta.tow.central %>% group_by(LIVE_CORAL) %>% count() %>% ungroup() %>% pull(LIVE_CORAL)
            lookup1 = lookup %>% filter(LIVE_CORAL %in% cats) 
            out2=sweep(pred1, 2, lookup1$Cover, '*')
            out3 = apply(out2, 3, rowSums)
            dat.central_brms.cumulative = ndata %>% cbind(tidyMCMC(as.mcmc(out3), conf.int=TRUE, conf.method='HPDinterval'))
            ## dat.central_brms.cumulative %>%
            ##     ggplot(aes(y=estimate, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.3) +
            ##     geom_line() +
            ##     geom_ribbon(data=dat.central, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.central, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.central.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.central.original, aes(y=mean), color='red') +
            ##     geom_line(data=b, aes(y=Mean, x=as.numeric(as.character(Year))), color='purple', size=2)
            save(mod.central_brms.cumulative, dat.central_brms.cumulative, file='../data/modelled/mod.central_brms.cumulative.RData')
        }  
    }
    ## ----end
    ## ---- Central.MGCV.tow.ordinal
    {
        if ('MGCV ordinal' %in% models) {
            library(mgcv)
            manta.tow.central = manta.tow %>%
                filter(Region=='Central GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL), nREEF_NAME=as.numeric(as.factor(REEF_NAME)))
            mod.central_mgcv.ordinal <- gam(nLIVE_CORAL ~ Year+s(nREEF_NAME, bs='re'),
                                            data=manta.tow.central,
                                            family=ocat(R=16),
                                            method='REML')

            lookup <- manta.tow.central %>% tidyr::expand(LIVE_CORAL) %>%
                mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))

            dat.central_mgcv.ordinal <- predict(mod.central_mgcv.ordinal, data.frame(Year=levels(manta.tow.central$Year), nREEF_NAME=0), type='response', exclude='s(nREEF_NAME)') %>% as.data.frame %>%
                setNames(c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>% 
                mutate(Year=levels(manta.tow.central$Year)) %>%
                pivot_longer(cols=-Year) %>%
                mutate(LIVE_CORAL=factor(name, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                         ordered=TRUE)) %>% 
                full_join(lookup) %>%
                mutate(P=value*Cover) %>%
                group_by(Year) %>%
                summarise(Mean=sum(P))

            ## dat.central.mgcv %>%
            ##     ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            ##     ## geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), alpha=0.3) +
            ##     geom_line() +  
            ##     geom_ribbon(data=dat.central, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.central, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.central.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.central.original, aes(y=mean), color='red')
            save(mod.central_mgcv.ordinal,dat.central_mgcv.ordinal, file='../data/modelled/mod.central_mgcv.ordinal.RData')
        }
    }
    ## ----end
    ## ---- Central.CLMM.tow.ordinal
    {
        if ('CLMM' %in% models) {
            library(ordinal) 
            manta.tow.central = manta.tow %>%
                filter(Region=='Central GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE))
            mod.central_clmm <-clmm(oLIVE_CORAL ~ Year + (1|REEF_NAME),
                                    data=manta.tow.central)
            
            summary(mod.central_clmm)
            dat.central_clmm=emmeans(mod.central_clmm, ~oLIVE_CORAL|Year, mode="prob") %>%
                as.data.frame() %>%
                group_by(Year) %>%
                mutate(oLIVE_CORAL=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>%
                ungroup %>%
                mutate(Perc=CoralTrends_calcPercent(oLIVE_CORAL),
                       Cover=prob*Perc) %>%
                group_by(Year) %>%
                summarise(Mean=sum(Cover))

            ## dat.central_clmm %>% 
            ##     ggplot(aes(x=as.numeric(as.character(Year)))) +
            ##     geom_line(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(data=dat.central, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.central, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.central.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.central.original, aes(y=mean), color='red')

            save(mod.central_clmm, dat.central_clmm, file='../data/modelled/dat.central_clmm.RData')
        }
    }
    ## ----end
    ## ---- Central.glmmTMB.tow.beta.linear
    {
        if ('glmmTMB.beta.linear' %in% models) {
            mod.central_glmmTMB.beta.linear <- glmmTMB(Cover ~ REPORT_YEAR + (1|REEF_NAME/REEF_YEAR),
                                                       data=manta.tow.central,
                                                       ## weights=dat.all.central$Tows,
                                                       family=beta_family())
            dat.central_glmmTMB.beta.linear = emmeans(mod.central_glmmTMB.beta.linear, ~REPORT_YEAR, at=list(REPORT_YEAR=unique(manta.tow.central$REPORT_YEAR)), type='response') %>%
                as.data.frame()

            save(mod.central_glmmTMB.beta.linear, dat.central_glmmTMB.beta.linear, file=paste0('../data/modelled/mod.central_glmmTMB.beta.linear.RData'))
            load(file=paste0('../data/modelled/mod.central_glmmTMB.beta.disp.RData'))
            dat.central_glmmTMB.beta.linear %>%
                ggplot() +
                geom_line(data=dat.central_glmmTMB.beta.disp, aes(y=response, x=as.numeric(as.character(Year))), color='blue') +
                geom_ribbon(data=dat.central_glmmTMB.beta.disp,aes(ymin=lower.CL, ymax=upper.CL, x=as.numeric(as.character(Year))), fill='lightblue', alpha=0.3) +
                geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, x=REPORT_YEAR), fill='orange', alpha=0.5) +
                geom_line(aes(y=response, x=REPORT_YEAR)) +
                scale_x_continuous('') +
                scale_y_continuous('Hard coral cover (%)', labels=function(x) x*100) +
                theme_classic() +
                ggtitle('Central GBR')
        }
    }
    ## ----end

    ## Compare models
    {
        if (COMPARE_MODELS) {        
            ## ---- beta.disp
        {
            load(file=paste0('../data/modelled/mod.central_glmmTMB.beta.disp.RData'))
            load(file=paste0('../data/modelled/mod.central_glmmTMB.beta.disp.rs.RData'))
            load(file=paste0('../data/modelled/mod.central_brms.beta.disp.RData'))
            load(file=paste0('../data/modelled/mod.central_inla.beta.disp.RData'))

            pdf(file = '../output/figures/comparison.central.pdf')
            g1 <- ggplot() +
                geom_line(data = dat.central_brms.beta.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'brms')) +
                geom_ribbon(data = dat.central_brms.beta.disp,
                            aes(y = response,
                                ymin = lower.HPD, ymax = upper.HPD,
                                x = as.numeric(as.character(Year)),
                                fill = 'brms'),
                            alpha=0.3) +
                geom_line(data = dat.central_glmmTMB.beta.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'glmmTMB')) +
                geom_ribbon(data = dat.central_glmmTMB.beta.disp,
                            aes(y = response,
                                ymin = lower.CL, ymax = upper.CL,
                                x = as.numeric(as.character(Year)),
                                fill = 'glmmTMB'),
                            alpha=0.3) + 
                geom_line(data = dat.central_inla.beta.disp,
                          aes(y = mean,
                              x = as.numeric(as.character(Year)),
                              colour = 'inla')) +
                geom_ribbon(data = dat.central_inla.beta.disp,
                            aes(y = mean,
                                ymin = lower, ymax = upper,
                                x = as.numeric(as.character(Year)),
                                fill = 'inla'),
                            alpha=0.3) 
            print(g1)
            dev.off()
        }
            ## ----end
            ## ---- beta.ry.disp
        {
            load(file=paste0('../data/modelled/mod.central_glmmTMB.beta.ry.disp.RData'))
            load(file=paste0('../data/modelled/mod.central_inla.beta.ry.disp.RData'))
            pdf(file = '../output/figures/comparison.central.ry.pdf')
            g1 <- ggplot() +
                ## geom_line(data = dat.central_brms.beta.ry.disp,
                ##           aes(y = response,
                ##               x = as.numeric(as.character(Year)),
                ##               colour = 'brms')) +
                ## geom_ribbon(data = dat.central_brms.beta.ry.disp,
                ##             aes(y = response,
                ##                 ymin = lower.HPD, ymax = upper.HPD,
                ##                 x = as.numeric(as.character(Year)),
                ##                 fill = 'brms'),
                ##             alpha=0.3) +
                geom_line(data = dat.central_glmmTMB.beta.ry.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'glmmTMB')) +
                geom_ribbon(data = dat.central_glmmTMB.beta.ry.disp,
                            aes(y = response,
                                ymin = lower.CL, ymax = upper.CL,
                                x = as.numeric(as.character(Year)),
                                fill = 'glmmTMB'),
                            alpha=0.3) + 
                geom_line(data = dat.central_inla.beta.ry.disp,
                          aes(y = mean,
                              x = as.numeric(as.character(Year)),
                              colour = 'inla')) +
                geom_ribbon(data = dat.central_inla.beta.ry.disp,
                            aes(y = mean,
                                ymin = lower, ymax = upper,
                                x = as.numeric(as.character(Year)),
                                fill = 'inla'),
                            alpha=0.3) 
            print(g1)
            dev.off()
        }
            ## ----end

        }

    }

    if (1==2) {
        ## Compare the models
    {
        ## ---- Central Compare models


        load(file='../data/modelled/dat.central.RData')
        dat.central.original <- dat.central
        load(file='../data/modelled/mod.central.RData')
        load(file='../data/modelled/mod.central_inla_beta.RData')
        load(file='../data/modelled/newdata.central_inla_beta.RData')
        newdata.central_beta <- newdata.central
        load(file='../data/modelled/mod.central_inla_binomial.RData')
        load(file='../data/modelled/newdata.central_inla_binomial.RData')
        newdata.central_binomial <- newdata.central
        load(file='../data/modelled/dat.central_inla_tow.RData')
        load(file='../data/modelled/dat.central_mgcv.RData')
        load(file='../data/modelled/dat.central_brms.cumulative.RData')
        load(file='../data/modelled/dat.central_clmm.RData')
        load(file='../data/modelled/dat.central_brms.beta.RData')
        load(file='../data/modelled/dat.central_glmmTMB.RData')
        load(file='../data/modelled/dat.central_mgcv.beta.RData')

        ## original
        g1 <- dat.central.original %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('Original (stan binomial reef level)') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g1

        ## inla (reef level binomial)
        g2 <- newdata.central_binomial %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA reef level binomial') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g2

        ## inla (reef level beta)
        g3 <- newdata.central_beta %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA reef level beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g3

        ## inla (tow level beta)
        g4 <- dat.central %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA tow level beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g4

        ## mgcv ordinal
        g5 <- dat.central.mgcv %>%
            ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('mgcv ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))

        ## brms ordinal
        g6 <- dat.central_brms.cumulative %>%
            ggplot(aes(y=estimate, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill='blue', alpha=0.3) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g6    

        ## clmm ordinal
        g7 <- dat.central_clmm %>%
            ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('clmm ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g7    

        ## brms beta
        g8 <- dat.central_brms.beta %>%
            ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g8

        library(patchwork)
        g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8


        dat.central_brms.beta %>%
            ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ## geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill='BRMS beta'), alpha=0.3) +
            geom_line(aes(color='BRMS beta')) +
            ## geom_ribbon(data=dat.central_brms.cumulative, aes(y=estimate, ymin=conf.low, ymax=conf.high, fill='BRMS ordinal'), alpha=0.3) +
            geom_line(data=dat.central_brms.cumulative, aes(y=estimate, color='BRMS ordinal')) +
            ## geom_ribbon(data=newdata.central_beta, aes(y=mean, ymin=lower, ymax=upper, fill='INLA beta'), alpha=0.3) +
            geom_line(data=newdata.central_beta, aes(y=mean, color='INLA beta')) +
            geom_line(data=dat.central_clmm, aes(y=Mean, color='clmm ordinal')) +
            geom_line(data=dat.central.mgcv, aes(y=Mean, color='mgcv ordinal')) +
            geom_line(data=dat.central_glmmTMB, aes(y=response, color='glmmTMB beta')) +
            geom_line(data=dat.central_mgcv.beta, aes(y=response, color='mgcv beta')) +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_brewer('Model', type='qual') + 
            scale_fill_discrete('Model') + 
            ## rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))


        rm(list=c('dat.central','mod.central','mod.central_inla_beta','newdata.central','newdata.central_beta', 'mod_central_inla_binomial','newdata.central_binomial'))
    }
        ## ----end
        ## ---- junk
    {
        ## summary(mod.central)
        ## dat.central = data.frame(Location='Central',Year=unique(dat.all.central$Year), N=length(unique(dat.all.central$REEF_NAME)))
        ## Xmat = model.matrix(~Year, dat.central)
        ## coefs = data.frame(mod.central) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
        ## coefs = mod.central$fit %>% as.data.frame() %>% dplyr::select(matches('b_.*'))
        ## Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
        ## dat.central = cbind(dat.central,
        ##             plyr:::adply(Fit,2,function(x) {
        ##     data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
        ##             })
        ## )


        ##     ## glmmTMB
        ##     library(glmmTMB)
        ##     library(emmeans)
        ##     mod.glmmTMB <-  glmmTMB(Cover ~ Year+(Year|REEF_NAME),
        ##                        data=dat.all.central,
        ##                        #weights=dat.all.gbr$Weight,
        ##                        family=beta_family)
        ##     emmeans(mod.glmmTMB, ~Year, type='response') %>% as.data.frame %>%
        ##     ggplot() +
        ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, x=as.numeric(as.character(Year))), fill='blue', alpha=0.2) +
        ##     geom_line(aes(y=response, x=as.numeric(as.character(Year)))) +
        ##     theme_bw()


        ## dat.central %>% arrange(Year)
        ## dat.central %>%
        ##     ggplot(aes(y=mean, x=as.numeric(Year))) +
        ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        ##     geom_line()

        ## save(dat.central, file='../data/modelled/dat.central.RData')
        ## save(dat.all.central, file='../data/modelled/dat.all.central.RData')
        ## save(mod.central, file='../data/modelled/mod.central.RData')
        ## rm(list=c('last_year','mod.central', 'dat.central', 'coefs', 'Fit'))
        ## gc()
    }
        ## ----end
    }
## ----end
}

## Southern
{
    ## ---- Southern
    ## ---- Southern.Data
    {
        dat.all.southern = dat.all %>%
            filter(Location=='Southern GBR') %>%
            droplevels %>% 
            mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
            group_by(REEF_NAME) %>%
            mutate(W=mean(Tows, na.rm=TRUE)) %>%
            ungroup %>%
            mutate(W1=W/sum(W)) %>%
            group_by(Year) %>%
            mutate(W2=Tows/sum(Tows)) %>%
            ungroup
        save(dat.all.southern, file='../data/modelled/dat.all.southern.RData') 

        ## Tow level data
        manta.tow.southern = manta.tow %>%
            filter(Region=='Southern GBR') %>%
            droplevels %>%
            mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                                      levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                      ordered=TRUE),
                   nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
                   nLIVE_CORAL=as.numeric(oLIVE_CORAL),
                   REEF_YEAR = interaction(REEF_NAME, Year)
                   )
        save(manta.tow.southern, file='../data/modelled/manta.tow.southern.RData')
    }
    ## ----end
    ## ---- Southern.Raw cells
    {
        dat.all.southern.cellmeans <- cellMeansRaw(dat.all.southern)
        rawAdd <- ggproto_Raw(dat.all.southern.cellmeans)

        manta.tow.southern.cellmeans <- cellMeansRaw(manta.tow %>%
                                                     filter(Region=='Northern GBR') %>%
                                                     droplevels %>% 
                                                     group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                                     summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                                     ungroup)
    }
    ## ----end                                        

    ## ---- Southern.stan_glmer.reef.binomial
    {
        if ('original' %in% models) {
            ## mod.southern <- ModelOriginal(form=formula(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME)),
            ##                               dat=dat.all.southern,
            ##                               location='Southern')
            ## mod.southern[['newdata']] %>%
            ##     ggplot(aes(y=mean, x=as.numeric(Year))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line()
            ## mod.southern[['mod']] %>% save(file='../data/modelled/mod.southern_original.RData')
            ## mod.southern[['newdata']] %>% save(file='../data/modelled/newdata.southern_original.RData')

            mod.southern_glmer <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                              data=dat.all.southern, family=binomial,iter=5000,warmup=2500,
                                              chains=3,cores=3)
            dat.southern_glmer <- data.frame(Location='Southern GBR',Year=unique(dat.all.southern$Year), N=length(unique(dat.all.southern$REEF_NAME)))
            Xmat <- model.matrix(~Year, dat.southern_glmer)

            coefs = data.frame(mod.southern_glmer) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
            Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
            dat.southern_glmer = cbind(dat.southern_glmer,
                                       plyr:::adply(Fit,2,function(x) {
                                           data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                                       })
                                       )
            save(dat.southern_glmer, mod.southern_glmer, file='../data/modelled/dat.southern_glmer.RData')
            ## save(dat.southern, file='../data/modelled/dat.southern.RData')
            ## save(mod.southern, file='../data/modelled/mod.southern.RData')
            rm(list=c('last_year','dat.southern_glmer','mod.southern_glmer'))
            gc()
            ## mod.southern <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
            ##                        data=dat.all.southern, family=binomial,iter=5000,warmup=2500,
            ##                        chains=3,cores=3)
            ## dat.southern <- data.frame(Location='Southern GBR',Year=unique(dat.all.southern$Year), N=length(unique(dat.all.southern$REEF_NAME)))
            ## Xmat <- model.matrix(~Year, dat.southern)

            ## coefs = data.frame(mod.southern) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
            ## Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
            ## dat.southern = cbind(dat.southern,
            ##                 plyr:::adply(Fit,2,function(x) {
            ##                     data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
            ##                 })
            ##                 )
            ## save(dat.southern, file='../data/modelled/dat.southern.RData')
            ## save(mod.southern, file='../data/modelled/mod.southern.RData')
            ## rm(list=c('last_year','mod.southern'))
            ## gc()
            ## mod.southern <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME),
            ##                             data=dat.all.southern,
            ##                             family=binomial,
            ##                             iter=5000,
            ##                             warmup=2500,
            ##                             chains=3,cores=3)
        }
    }
    ## ----end
    ## ---- Southern.stan_glmer.reef.beta
    {
        if ('stan_glmer beta' %in% models) {
            mod.southern <- Modelstan_glmer(form=Cover ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                            dat=dat.all.southern)
            newdata.southern <- cellMeansOriginal(mod.southern, dat=dat.all.southern, location='Southern') 
            newdata.southern %>%
                ggplot(aes(y=mean, x=as.numeric(Year))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                geom_line()
            save(mod.southern, file='../data/modelled/mod.southern_stan_glmer.RData')
            save(newdata.southern, file='../data/modelled/newdata.southern_stan_glmer.RData')
            rm(list=c('last_year','mod.southern', 'newdata.southern', 'coefs', 'Fit'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.BRMS.reef.beta
    {
        if ('BRMS_reef beta' %in% models) {
            priors <- prior(normal(0, 5), class = "b") +
                prior(normal(0, 5), class = "Intercept") +
                ## prior(gamma(2, 0.1), class = "sd") +
                prior(gamma(1, 0.5), class = "sd") +
                prior(gamma(0.01, 0.01), class = "phi")
            inits = list(list(phi=list(rgamma(1,0.1,0.1))),
                         list(phi=list(rgamma(1,0.1,0.1))),
                         list(phi=list(rgamma(1,0.1,0.1)))
                         )
            mod.southern <- Modelbrms(form=Cover|weights(Tows) ~ Year + (1|REEF_NAME),
                                      dat=dat.all.southern)
            
            mod.southern <- Modelbrms(form=Cover|weights(Tows) ~ Year + (Year|REEF_NAME),
                                      dat=dat.all.southern)
        }
    }
    ## ----end
    ## ---- Southern.INLA.reef.beta vanilla
    {
        if ("INLA_reef beta" %in% models) {
            dat.inla <- dataINLA(dat=dat.all.southern %>% mutate(W=NA))
            dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
            mod.southern_inla_reef.beta <- ModelINLA_beta(form=Cover~Year+
                                                              f(P_CODE.mod, model='iid')+
                                                              f(REEF_NAME, model='iid'),
                                                          dat=dat.inla[['dat.1']])
            dat.southern_inla_reef.beta <- cellMeansINLA(mod=mod.southern_inla_reef.beta, newdata.hcc=dat.inla[['newdata.hcc']],
                                                         n.2=dat.inla[['n.2']])
            ## newdata.southern %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd
            save(dat.southern_inla_reef.beta, mod.southern_inla_reef.beta, file='../data/modelled/mod.southern_inla_reef.beta.RData')
            ## save(newdata.southern, file='../data/modelled/newdata.southern_inla_beta.RData')
            rm(list=c('mod.southern_inla_reef.beta', 'dat.southern_inla_reef.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.INLA.reef.beta scaled
    {
        if ("INLA_reef beta" %in% models) {
            dat.inla <- dataINLA(dat=dat.all.southern %>% mutate(W=NA))
            dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
            mod.southern_inla_reef.beta.scaled <- ModelINLA_beta(form=Cover~Year+
                                                                     f(P_CODE.mod, model='iid')+
                                                                     f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                                                     f(REEF_NAME1, Year, model='iid', scale=dat.scale$Tows),
                                                                 dat=dat.inla[['dat.1']])
            dat.southern_inla_reef.beta.scaled <- cellMeansINLA(mod=mod.southern_inla_reef.beta.scaled, newdata.hcc=dat.inla[['newdata.hcc']],
                                                                n.2=dat.inla[['n.2']])
            ## newdata.southern %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd
            save(dat.southern_inla_reef.beta.scaled, mod.southern_inla_reef.beta.scaled, file='../data/modelled/mod.southern_inla_reef.beta.scaled.RData')
            ## save(newdata.southern, file='../data/modelled/newdata.southern_inla_beta.RData')
            rm(list=c('mod.southern_inla_reef.beta.scaled', 'dat.southern_inla_reef.beta.scaled'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.INLA.reef.binomial
    {
        if ('INLA_reef binomial' %in% models) {
            dat.inla <- dataINLA(dat=dat.all.southern)
            mod.southern <- ModelINLA_binomial(form=Cvr1~Year+
                                                   f(P_CODE.mod, model='iid')+
                                                   f(REEF_NAME, model='iid') +
                                                   f(REEF_NAME1, Year, model='iid'),
                                               dat=dat.inla[['dat.1']])
            newdata.southern <- cellMeansINLA(mod=mod.southern, newdata.hcc=dat.inla[['newdata.hcc']],
                                              n.2=dat.inla[['n.2']])
            newdata.southern %>%
                ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
                geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
                geom_line() +
                rawAdd
            save(mod.southern, file='../data/modelled/mod.southern_inla_binomial.RData')
            save(newdata.southern,file='../data/modelled/newdata.southern_inla_binomial.RData')
            rm(list=c('mod.southern', 'newdata.southern'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.glmmTMB.reef.beta
    {
        if ('glmmTMB_reef beta' %in% models) {
            mod.southern_glmmTMB.reef.beta <- glmmTMB(Cover ~ Year + (1|P_CODE.mod) +
                                                          (1|REEF_NAME),
                                                      data=dat.all.southern,
                                                      weights=dat.all.southern$Tows,
                                                      family=beta_family())
            dat.southern_glmmTMB.reef.beta  = emmeans(mod.southern_glmmTMB.reef.beta, ~Year, type='response') %>%
                as.data.frame()
            save(mod.southern_glmmTMB.reef.beta, dat.southern_glmmTMB.reef.beta, file='../data/modelled/dat.southern_glmmTMB.reef.beta.RData')
            
            ## emmeans(mod.southern, ~Year, type='response') %>%
            ##     as.data.frame() %>% 
            ##     ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd

            rm(list=c('dat.southern_glmmTMB.reef.beta', 'mod.southern_glmmTMB.reef.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.glmmTMB.reef.beta
    {
        if ('glmmTMB_reef beta' %in% models) {
            mod.southern_glmmTMB.reef.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME),
                                                           dispformula = ~Year,
                                                           data=dat.all.southern,
                                                           weights=dat.all.southern$Tows,
                                                           family=beta_family())
            dat.southern_glmmTMB.reef.beta.disp  = emmeans(mod.southern_glmmTMB.reef.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            save(mod.southern_glmmTMB.reef.beta.disp, dat.southern_glmmTMB.reef.beta.disp, file='../data/modelled/dat.southern_glmmTMB.reef.beta.disp.RData')
            
            ## emmeans(mod.southern, ~Year, type='response') %>%
            ##     as.data.frame() %>% 
            ##     ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
            ##     geom_line() +
            ##     rawAdd

            rm(list=c('dat.southern_glmmTMB.reef.beta.disp', 'mod.southern_glmmTMB.reef.beta.disp'))
            gc()
        }
    }
    ## ----end

    ## ---- Southern.INLA.tow.beta scaled
    {
        if ("INLA_tow beta scaled" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.southern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.southern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            mod.southern_inla.beta.scaled <- ModelINLA_beta(form=Cover~Year +
                                                                f(P_CODE.mod, model='iid') +
                                                                f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                                                f(REEF_NAME2, YEAR1, model='iid', scale=dat.scale$Tows),
                                                            dat=dat,
                                                            family='beta',
                                                            weights=NULL)
            
            dat.southern_inla.beta.scaled <- cellMeansINLA(mod=mod.southern_inla.beta.scaled, newdata.hcc=dat.inla[['newdata.hcc']],
                                                           n.2=dat.inla[['n.2']], FUN=plogis)
            ## dat.raw = dat %>% group_by(REEF_NAME, Year) %>% summarize(Cover=mean(Cover, na.rm=TRUE)) %>% filter(!is.na(Cover)) %>% droplevels
            ## dat.southern %>%
            ##     ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill='blue') +
            ##     geom_line(color='blue') +
            ##    rawAdd +
            ##    ggtitle('GBR INLA tow level') #
            ##    ## geom_line(data=mod.brms.df, aes(y=response, x=as.numeric(as.character(year)))) +
            ##    ## geom_ribbon(data=mod.brms.df, aes(y=response, ymin=lower.HPD, ymax=upper.HPD, x=as.numeric(as.character(year))), alpha=0.3)        

            save(mod.southern_inla.beta.scaled, dat.southern_inla.beta.scaled, file='../data/modelled/mod.southern_inla.beta.scaled.RData')
            rm(list=c('dat.southern_inla.beta.scaled', 'mod.southern_inla.beta.scaled'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.INLA.tow.beta
    {
        if ("INLA_tow beta" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.southern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.southern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            mod.southern_inla.beta <- ModelINLA_beta(form=Cover~Year +
                                                         f(P_CODE.mod, model='iid') +
                                                         f(REEF_NAME, model='iid') +
                                                         f(REEF_YEAR, model='iid'), 
                                                     dat=dat,
                                                     family='beta',
                                                     weights=NULL)
            
            dat.southern_inla.beta <- cellMeansINLA(mod=mod.southern_inla.beta, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.southern_inla.beta, dat.southern_inla.beta, file='../data/modelled/mod.southern_inla.beta.RData')
            rm(list=c('dat.southern_inla.beta', 'mod.southern_inla.beta'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.INLA.tow.beta **
    {
        if ("INLA_tow beta disp" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.southern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.southern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            ## we will leverage a mixed likelihood model
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            
            dd <- dat %>%
                dplyr::select(Cover, Year, REEF_NAME) %>%
                mutate(YEAR = Year)%>%
                pivot_wider(id_cols = c(YEAR,REEF_NAME),
                            names_from = Year,
                            values_from = Cover)
            dd1 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(matches("[0-9]{4}")) %>%
                as.matrix()
            dd2 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(YEAR, REEF_NAME)

            mod.southern_inla.beta.disp <- inla(form = dd1~YEAR +
                                               f(REEF_NAME, model='iid'),
                                           dat=dd2,
                                           family=rep('beta',ncol(dd1)),
                                           control.fixed = list(mean = 0, prec = 0.001,
                                                                mean.intercept = 0.5,
                                                                prec.intercept = 0.001),
                                           control.predictor = list(compute = TRUE,
                                                                    link = 1,
                                                                    quantiles = c(0.025,0.25,0.5,0.75,0.975)
                                                                    )
                                           )
            dat.southern_inla.beta.disp <- cellMeansINLA(mod=mod.southern_inla.beta.disp, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.southern_inla.beta.disp, dat.southern_inla.beta.disp, file='../data/modelled/mod.southern_inla.beta.disp.RData')
            rm(list=c('dat.southern_inla.beta.disp', 'mod.southern_inla.beta.disp'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.INLA.tow.ry.beta **
    {
        if ("INLA_tow beta ry disp" %in% models) {
            dat.inla <- dataINLA(dat=manta.tow.southern, level='tow')
            dat = dat.inla[['dat.1']]
            dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
            rawAdd <- ggproto_Raw(dat.cellmeans)
            dat.scale = manta.tow.southern %>% group_by(REEF_NAME, REPORT_YEAR) %>%
                summarise(Tows=length(TOW_SEQ_NO)) %>% ungroup %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows))
            ## we will leverage a mixed likelihood model
            dat = dat %>% mutate(REEF_YEAR = interaction(REEF_NAME, Year))
            
            dd <- dat %>%
                dplyr::select(Cover, Year, REEF_NAME, REEF_YEAR) %>%
                mutate(YEAR = Year)%>%
                pivot_wider(id_cols = c(YEAR,REEF_NAME, REEF_YEAR),
                            names_from = Year,
                            values_from = Cover)
            dd1 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(matches("[0-9]{4}")) %>%
                as.matrix()
            dd2 <- dd %>% unnest(matches("[0-9]{4}")) %>%
                dplyr::select(YEAR, REEF_NAME, REEF_YEAR)

            mod.southern_inla.beta.ry.disp <- inla(form = dd1~YEAR +
                                               f(REEF_NAME, model='iid') +
                                               f(REEF_YEAR, model='iid'),
                                           dat=dd2,
                                           family=rep('beta',ncol(dd1)),
                                           control.fixed = list(mean = 0, prec = 0.001,
                                                                mean.intercept = 0.5,
                                                                prec.intercept = 0.001),
                                           control.predictor = list(compute = TRUE,
                                                                    link = 1,
                                                                    quantiles = c(0.025,0.25,0.5,0.75,0.975)
                                                                    )
                                           )
            dat.southern_inla.beta.ry.disp <- cellMeansINLA(mod=mod.southern_inla.beta.ry.disp, newdata.hcc=dat.inla[['newdata.hcc']],
                                                    n.2=dat.inla[['n.2']], FUN=plogis)
            save(mod.southern_inla.beta.ry.disp, dat.southern_inla.beta.ry.disp, file='../data/modelled/mod.southern_inla.beta.ry.disp.RData')
            rm(list=c('dat.southern_inla.beta.ry.disp', 'mod.southern_inla.beta.ry.disp'))
            gc()
        }
    }
    ## ----end
    ## ---- Southern.BRMS.tow.beta vanilla
    {
        if ('BRMS beta vanilla' %in% models) {
            mod.southern_brms.beta <- brm(Cover ~ Year + (1|REEF_NAME),
                                          data=manta.tow.southern,
                                          family=Beta(link='logit'),
                                          iter=1e4,
                                          warmup=5e3,
                                          thin=5,
                                          chains=4, cores=4,
                                          prior = prior(normal(0, 3), class = "b") +
                                              prior(normal(0, 3), class = "Intercept") +
                                              prior(gamma(2, 1), class = "sd") +
                                              prior(gamma(2, 1), class = "sd", group = "REEF_NAME") +
                                              prior(gamma(2, 1), class = "phi")
                                          )
            dat.southern_brms.beta = emmeans(mod.southern_brms.beta, ~Year, type='response') %>%
                as.data.frame()
            save(mod.southern_brms.beta, dat.southern_brms.beta, file=paste0('../data/modelled/mod.southern_brms_beta.RData'))
            ## dat.southern_brms.beta %>%
            ##     ggplot(aes(y=response, x=as.numeric(as.character(year)))) +
            ##     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), alpha=0.3) +
            ##     geom_line() +
            ##     geom_ribbon(data=dat.southern, aes(y=mean, x=as.numeric(as.character(Year)), ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.southern, aes(y=mean, x=as.numeric(as.character(Year))), color='blue') +
            ##     geom_ribbon(data=dat.southern.original, aes(y=mean, x=as.numeric(as.character(Year)), ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.southern.original, aes(y=mean, x=as.numeric(as.character(Year))), color='red')
            ## save(dat.southern_brms.beta, file='../data/modelled/dat.southern_brms.beta.RData')
            rm(list=c('dat.southern_brms.beta', 'mod.southern_brms.beta'))
            gc()
        }    
    }
    ## ----end
    ## ---- Southern.BRMS.tow.beta disp **
    {
        if ('BRMS beta disp' %in% models & 'southern' %in% zone) {
            priors <- prior(normal(0, 3), class = "b") +
                prior(normal(0, 3), class = "Intercept") +
                prior(gamma(2, 1), class = "sd") 
            ## The above priors where 0,1  0,1  2,1
            ## might like to try 0,2 0,1.5 2,1
            mod.southern_brms.beta.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME), phi~0+Year),
                                               data=manta.tow.southern,
                                               family=Beta(link='logit'),
                                               iter=1e4,
                                               warmup=5e3,
                                               thin=5,
                                               chains=4, cores=4,
                                               prior = priors
                                               ## prior = prior(normal(0, 3), class = "b") +
                                               ##     prior(normal(0, 3), class = "Intercept") +
                                               ##     prior(gamma(2, 1), class = "sd") #+
                                               ## ## prior(gamma(2, 1), class = "phi")
                                               )
            dat.southern_brms.beta.disp = emmeans(mod.southern_brms.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Southern.BRMS.tow.beta disp diagnostics
            {
                ## sampling diagnostics
                pdf(file = '../output/figures/traceplots_southern_brms.beta.disp.pdf')
                rstan::traceplot(mod.southern_brms.beta.disp$fit)
                dev.off()

                ## density overlay
                pdf(file = '../output/figures/density_southern_brms.beta.disp.pdf')
                mod.southern_brms.beta.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100) 
                dev.off()

                ## DHARMa residuals
                preds <- mod.southern_brms.beta.disp %>%
                    posterior_predict(ndraws = 250, summary = FALSE)
                mod.resids <- createDHARMa(
                    simulatedResponse = t(preds),
                    observedResponse = manta.tow.southern$Cover,
                    fittedPredictedResponse = apply(preds, 2, median),
                    integerResponse = FALSE
                )
                pdf(file = '../output/figures/DHARMa_southern_brms.beta.disp.pdf')
                mod.resids %>% plot()
                dev.off()
                save(mod.resids, file=paste0('../data/modelled/resids.southern_brms.beta.disp.RData'))
            }
            ## ----end
            save(mod.southern_brms.beta.disp, dat.southern_brms.beta.disp, file=paste0('../data/modelled/mod.southern_brms.beta.disp.RData'))
            rm(list=c('dat.southern_brms.beta.disp', 'mod.southern_brms.beta.disp'))
            gc()
        }  
    }
    ## ----end
    ## ---- Southern.BRMS.tow.beta disp **
    {
        if ('BRMS beta ry disp' %in% models & 'southern' %in% zone) {
            priors <- prior(normal(0, 3), class = "b") +
                prior(normal(0, 3), class = "Intercept") +
                prior(gamma(2, 1), class = "sd") 
            ## The above priors where 0,1  0,1  2,1
            ## might like to try 0,2 0,1.5 2,1
            mod.southern_brms.beta.ry.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME/REEF_YEAR), phi~0+Year),
                                               data=manta.tow.southern,
                                               family=Beta(link='logit'),
                                               iter=1e4,
                                               warmup=5e3,
                                               thin=5,
                                               chains=4, cores=4,
                                               prior = priors
                                               ## prior = prior(normal(0, 3), class = "b") +
                                               ##     prior(normal(0, 3), class = "Intercept") +
                                               ##     prior(gamma(2, 1), class = "sd") #+
                                               ## ## prior(gamma(2, 1), class = "phi")
                                               )
            dat.southern_brms.beta.ry.disp = emmeans(mod.southern_brms.beta.ry.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Southern.BRMS.tow.beta ry disp diagnostics
            {
                ## sampling diagnostics
                pdf(file = '../output/figures/traceplots_southern_brms.beta.ry.disp.pdf')
                rstan::traceplot(mod.southern_brms.beta.ry.disp$fit)
                dev.off()

                ## density overlay
                pdf(file = '../output/figures/density_southern_brms.beta.ry.disp.pdf')
                mod.southern_brms.beta.ry.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100) 
                dev.off()

                ## DHARMa residuals
                preds <- mod.southern_brms.beta.ry.disp %>%
                    posterior_predict(ndraws = 250, summary = FALSE)
                mod.resids <- createDHARMa(
                    simulatedResponse = t(preds),
                    observedResponse = manta.tow.southern$Cover,
                    fittedPredictedResponse = apply(preds, 2, median),
                    integerResponse = FALSE
                )
                pdf(file = '../output/figures/DHARMa_southern_brms.beta.ry.disp.pdf')
                mod.resids %>% plot()
                dev.off()
                save(mod.resids, file=paste0('../data/modelled/resids.southern_brms.beta.ry.disp.RData'))
            }
            ## ----end
            save(mod.southern_brms.beta.ry.disp, dat.southern_brms.beta.ry.disp, file=paste0('../data/modelled/mod.southern_brms.beta.ry.disp.RData'))
            rm(list=c('dat.southern_brms.beta.ry.disp', 'mod.southern_brms.beta.ry.disp'))
            gc()
        }  
    }
    ## ----end
    ## ---- Southern.MGCV.tow.beta
    {
        if ('MGCV beta' %in% models) {
            manta.tow.southern = manta.tow %>%
                filter(Region=='Southern GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL), nREEF_NAME=as.numeric(as.factor(REEF_NAME)))
            mod.southern_mgcv.beta <- gam(Cover ~ Year+s(nREEF_NAME, bs='re'),
                                          data=manta.tow.southern,
                                          family=betar,
                                          method='REML')

            dat.southern_mgcv.beta <- emmeans(mod.southern_mgcv.beta, ~Year, type='response') %>% as.data.frame() 
            save(dat.southern_mgcv.beta, mod.southern_mgcv.beta, file='../data/modelled/mod.southern_mgcv.beta.RData')
        }
    }
    ## ----end
    ## ---- Southern.glmmTMB.tow.beta vanilla
    {
        if ('glmmTMB_tow beta vanilla' %in% models) {
            mod.southern_glmmTMB.beta <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                                 data=manta.tow.southern,
                                                 ## weights=dat.all.southern$Tows,
                                                 family=beta_family())
            dat.southern_glmmTMB.beta = emmeans(mod.southern_glmmTMB.beta, ~Year, type='response') %>%
                as.data.frame()
            ## DHARMa::simulateResiduals(mod.southern_glmmTMB.beta, plot=TRUE)
            ## performance::check_model(mod.southern_glmmTMB.beta)
            save(mod.southern_glmmTMB.beta, dat.southern_glmmTMB.beta, file=paste0('../data/modelled/mod.southern_glmmTMB.beta.RData'))
        }
    }
    ## ----end
    ## ---- Southern.glmmTMB.tow.beta disp **
    {
        if ('glmmTMB_tow beta disp' %in% models & 'southern' %in% zone) {
            cat('Fitting brms disp for Southern\n\n')
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.southern_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + #(1|REEF_NAME/REEF_YEAR),
                                                          (1|REEF_NAME),
                                                      dispformula = ~Year,
                                                      data=manta.tow.southern,
                                                      ## weights=dat.all.southern$Tows,
                                                      family=beta_family(),
                                                     control = glmmTMBControl(parallel = nt))
            dat.southern_glmmTMB.beta.disp = emmeans(mod.southern_glmmTMB.beta.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Southern.glmmTMB.tow.beta disp diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_southern_glmmTMB.beta.disp.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.southern_glmmTMB.beta.disp,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.southern_glmmTMB.beta.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.southern_glmmTMB.beta.disp.RData'))
            }
            ## ----end
            ## DHARMa::simulateResiduals(mod.southern_glmmTMB.beta.disp, plot=TRUE)
            ## performance::check_model(mod.southern_glmmTMB.beta.disp)
            save(mod.southern_glmmTMB.beta.disp, dat.southern_glmmTMB.beta.disp, file=paste0('../data/modelled/mod.southern_glmmTMB.beta.disp.RData'))
        }
    }
    ## ----end
    ## ---- Southern.glmmTMB.tow.beta ry disp **
    {
        if ('glmmTMB_tow beta ry disp' %in% models & 'southern' %in% zone) {
            cat('Fitting brms ry disp for Southern\n\n')
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.southern_glmmTMB.beta.ry.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                                      dispformula = ~Year,
                                                      data=manta.tow.southern,
                                                      ## weights=dat.all.southern$Tows,
                                                      family=beta_family(),
                                                     control = glmmTMBControl(parallel = nt))
            dat.southern_glmmTMB.beta.ry.disp = emmeans(mod.southern_glmmTMB.beta.ry.disp, ~Year, type='response') %>%
                as.data.frame()
            ## ---- Southern.glmmTMB.tow.beta disp diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_southern_glmmTMB.beta.ry.disp.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.southern_glmmTMB.beta.ry.disp,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.southern_glmmTMB.beta.ry.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.southern_glmmTMB.beta.ry.disp.RData'))
            }
            ## ----end
            ## DHARMa::simulateResiduals(mod.southern_glmmTMB.beta.ry.disp, plot=TRUE)
            ## performance::check_model(mod.southern_glmmTMB.beta.ry.disp)
            save(mod.southern_glmmTMB.beta.ry.disp, dat.southern_glmmTMB.beta.ry.disp, file=paste0('../data/modelled/mod.southern_glmmTMB.beta.ry.disp.RData'))
        }
    }
    ## ----end
    ## ---- Southern.glmmTMB.tow.beta disp random.effects
    {
        if ('glmmTMB_tow beta disp' %in% models) {
            nt <- parallel::detectCores()
            nt <- ifelse(nt>5, 20, nt)
            mod.southern_glmmTMB.beta.disp.rs <- glmmTMB(Cover ~ Year + (Year|REEF_NAME),
                                                         dispformula = ~Year,
                                                         data=manta.tow.southern,
                                                         ## weights=dat.all.southern$Tows,
                                                         family=beta_family(),
                                                         control = glmmTMBControl(parallel=nt))
            dat.southern_glmmTMB.beta.disp.rs = emmeans(mod.southern_glmmTMB.beta.disp.rs, ~Year, type='response') %>%
                as.data.frame()
            ## did not converge
            ## ---- Southern.glmmTMB.tow.beta disp random.effects diagnostics
            {
                ## DHARMa residuals
                pdf(file = '../output/figures/DHARMa_southern_glmmTMB.beta.disp.rs.pdf')
                mod.resids <- DHARMa::simulateResiduals(mod.southern_glmmTMB.beta.disp.rs,
                                                       plot=TRUE, integer = TRUE)
                mod.resids %>% plot()
                dev.off()
                ## performance::check_model(mod.southern_glmmTMB.beta.disp)
                save(mod.resids, file=paste0('../data/modelled/resids.southern_glmmTMB.beta.disp.RData'))
                ## performance::check_model(mod.southern_glmmTMB.beta.disp.rs)
            }
            ## ----end
            save(mod.southern_glmmTMB.beta.disp.rs, dat.southern_glmmTMB.beta.disp.rs, file=paste0('../data/modelled/mod.southern_glmmTMB.beta.disp.rs.RData'))
        }
    }
    ## ----end
    ## ---- Southern.BRMS.tow.ordinal
    {
        if ('BRMS ordinal' %in% models) {
            library(brms)
            manta.tow.southern.tally <- manta.tow.southern %>%
                group_by(REEF_NAME, Year, oLIVE_CORAL) %>%
                count() %>%
                ungroup
            mod.southern_brms.cumulative <- brm(bf(oLIVE_CORAL|weights(n) ~ Year+(1|REEF_NAME)),
                                                data=manta.tow.southern.tally,
                                                ## family=cumulative("logit", threshold = 'flexible'),
                                                family=cumulative("logit", threshold = 'equidistant'),
                                                ## family=cumulative("probit", threshold = 'equidistant'),
                                                iter=1e4,
                                                warmup=5e3,
                                                thin=5,
                                                chains=4, cores=4,
                                                prior = prior(normal(0, 3), class = "b") +
                                                    prior(gamma(1, 0.5), class = "sd") +
                                                    prior(normal(0, 3), class = "Intercept") 
                                                )
            summary(mod.southern_brms.cumulative)
            ## save(mod.southern_brms.cumulative, file='../data/modelled/mod.southern_brms.cumulative.RData')

            ndata = manta.tow.southern %>% tidyr::expand(Year=levels(Year))
            pred1 <- posterior_epred(mod.southern_brms.cumulative, newdata=ndata, re_formula=NA) %>%
                aperm(c(1,3,2))
            lookup <-data.frame(LIVE_CORAL=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>%
                mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))
            ## lookup <- manta.tow.southern %>% tidyr::expand(LIVE_CORAL) %>%
            ##     mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))
            cats <- manta.tow.southern %>% group_by(LIVE_CORAL) %>% count() %>% ungroup() %>% pull(LIVE_CORAL)
            lookup1 = lookup %>% filter(LIVE_CORAL %in% cats) 
            out2=sweep(pred1, 2, lookup1$Cover, '*')
            out3 = apply(out2, 3, rowSums)
            dat.southern_brms.cumulative = ndata %>% cbind(tidyMCMC(as.mcmc(out3), conf.int=TRUE, conf.method='HPDinterval'))
            ## dat.southern_brms.cumulative %>%
            ##     ggplot(aes(y=estimate, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.3) +
            ##     geom_line() +
            ##     geom_ribbon(data=dat.southern, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.southern, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.southern.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.southern.original, aes(y=mean), color='red') +
            ##     geom_line(data=b, aes(y=Mean, x=as.numeric(as.character(Year))), color='purple', size=2)
            save(mod.southern_brms.cumulative, dat.southern_brms.cumulative, file='../data/modelled/mod.southern_brms.cumulative.RData')
        }  
    }
    ## ----end
    ## ---- Southern.MGCV.tow.ordinal
    {
        if ('MGCV ordinal' %in% models) {
            library(mgcv)
            manta.tow.southern = manta.tow %>%
                filter(Region=='Southern GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE),
                       nLIVE_CORAL=as.numeric(oLIVE_CORAL), nREEF_NAME=as.numeric(as.factor(REEF_NAME)))
            mod.southern_mgcv.ordinal <- gam(nLIVE_CORAL ~ Year+s(nREEF_NAME, bs='re'),
                                             data=manta.tow.southern,
                                             family=ocat(R=16),
                                             method='REML')

            lookup <- manta.tow.southern %>% tidyr::expand(LIVE_CORAL) %>%
                mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL))

            dat.southern_mgcv.ordinal <- predict(mod.southern_mgcv.ordinal, data.frame(Year=levels(manta.tow.southern$Year), nREEF_NAME=0), type='response', exclude='s(nREEF_NAME)') %>% as.data.frame %>%
                setNames(c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>% 
                mutate(Year=levels(manta.tow.southern$Year)) %>%
                pivot_longer(cols=-Year) %>%
                mutate(LIVE_CORAL=factor(name, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                                         ordered=TRUE)) %>% 
                full_join(lookup) %>%
                mutate(P=value*Cover) %>%
                group_by(Year) %>%
                summarise(Mean=sum(P))

            ## dat.southern.mgcv %>%
            ##     ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            ##     ## geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), alpha=0.3) +
            ##     geom_line() +  
            ##     geom_ribbon(data=dat.southern, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.southern, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.southern.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.southern.original, aes(y=mean), color='red')
            save(mod.southern_mgcv.ordinal,dat.southern_mgcv.ordinal, file='../data/modelled/mod.southern_mgcv.ordinal.RData')
        }
    }
    ## ----end
    ## ---- Southern.CLMM.tow.ordinal
    {
        if ('CLMM' %in% models) {
            library(ordinal) 
            manta.tow.southern = manta.tow %>%
                filter(Region=='Southern GBR') %>%
                droplevels %>%
                mutate(oLIVE_CORAL=factor(LIVE_CORAL, levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'), ordered=TRUE))
            mod.southern_clmm <-clmm(oLIVE_CORAL ~ Year + (1|REEF_NAME),
                                     data=manta.tow.southern)
            
            summary(mod.southern_clmm)
            dat.southern_clmm=emmeans(mod.southern_clmm, ~oLIVE_CORAL|Year, mode="prob") %>%
                as.data.frame() %>%
                group_by(Year) %>%
                mutate(oLIVE_CORAL=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U')) %>%
                ungroup %>%
                mutate(Perc=CoralTrends_calcPercent(oLIVE_CORAL),
                       Cover=prob*Perc) %>%
                group_by(Year) %>%
                summarise(Mean=sum(Cover))

            ## dat.southern_clmm %>% 
            ##     ggplot(aes(x=as.numeric(as.character(Year)))) +
            ##     geom_line(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            ##     geom_ribbon(data=dat.southern, aes(y=mean, ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            ##     geom_line(data=dat.southern, aes(y=mean), color='blue') +
            ##     geom_ribbon(data=dat.southern.original, aes(y=mean, ymin=lower, ymax=upper), fill='red', alpha=0.3) +
            ##     geom_line(data=dat.southern.original, aes(y=mean), color='red')

            save(mod.southern_clmm, dat.southern_clmm, file='../data/modelled/dat.southern_clmm.RData')
        }
    }
    ## ----end
    ## ---- Southern.glmmTMB.tow.beta.linear
    {
        if ('glmmTMB.beta.linear' %in% models) {
            mod.southern_glmmTMB.beta.linear <- glmmTMB(Cover ~ REPORT_YEAR + (1|REEF_NAME/REEF_YEAR),
                                                        data=manta.tow.southern,
                                                        ## weights=dat.all.southern$Tows,
                                                        family=beta_family())
            dat.southern_glmmTMB.beta.linear = emmeans(mod.southern_glmmTMB.beta.linear, ~REPORT_YEAR, at=list(REPORT_YEAR=unique(manta.tow.southern$REPORT_YEAR)), type='response') %>%
                as.data.frame()

            save(mod.southern_glmmTMB.beta.linear, dat.southern_glmmTMB.beta.linear, file=paste0('../data/modelled/mod.southern_glmmTMB.beta.linear.RData'))
            load(file=paste0('../data/modelled/mod.southern_glmmTMB.beta.disp.RData'))
            dat.southern_glmmTMB.beta.linear %>%
                ggplot() +
                geom_line(data=dat.southern_glmmTMB.beta.disp, aes(y=response, x=as.numeric(as.character(Year))), color='blue') +
                geom_ribbon(data=dat.southern_glmmTMB.beta.disp,aes(ymin=lower.CL, ymax=upper.CL, x=as.numeric(as.character(Year))), fill='lightblue', alpha=0.3) +
                geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, x=REPORT_YEAR), fill='orange', alpha=0.5) +
                geom_line(aes(y=response, x=REPORT_YEAR)) +
                scale_x_continuous('') +
                scale_y_continuous('Hard coral cover (%)', labels=function(x) x*100) +
                theme_classic() +
                ggtitle('Southern GBR')

            emmeans(mod.southern_glmmTMB.beta.linear, ~REPORT_YEAR, at=list(REPORT_YEAR=c(2000:2002, 2010:2011))) %>%
                regrid() %>%
                confint() %>%
                as.data.frame()
            ## absolute annual change
            ## coral cover declines by 0.345 (e.g from 30.446% to 30.102%) 
            emmeans(mod.southern_glmmTMB.beta.linear, ~REPORT_YEAR, at=list(REPORT_YEAR=2010:2012)) %>%
                regrid() %>%
                pairs(reverse=TRUE) %>%
                summary(infer=c(TRUE,TRUE)) %>%
                dplyr::mutate(across(c('estimate',ends_with('.CL')), function(x) x*100))
            ## confint()
            ## percentage annual change
            ## coral cover declines by 1.33% per year.
            ## this could be as great as a 1.47% decline or as mild as a 0.8% decline
            emmeans(mod.southern_glmmTMB.beta.linear, ~REPORT_YEAR, at=list(REPORT_YEAR=2000:2002)) %>%
                regrid(transform='log') %>%
                pairs(reverse = TRUE) %>%
                summary(infer=c(TRUE,TRUE)) %>%
                dplyr::mutate(across(c('estimate',ends_with('.CL')), function(x) -100*(1-exp(x))))
            ## confint() %>%
            ## dplyr::mutate(across(where(is.numeric), exp))
        }
    }
    ## ----end

    ## Compare models
    {
        if (COMPARE_MODELS) {        
            ## ---- beta.disp
            {
            load(file=paste0('../data/modelled/mod.southern_glmmTMB.beta.disp.RData'))
            load(file=paste0('../data/modelled/mod.southern_glmmTMB.beta.disp.rs.RData'))
            load(file=paste0('../data/modelled/mod.southern_brms.beta.disp.RData'))
            load(file=paste0('../data/modelled/mod.southern_inla.beta.disp.RData'))

            pdf(file = '../output/figures/comparison.southern.pdf')
            g1 <- ggplot() +
                geom_line(data = dat.southern_brms.beta.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'brms')) +
                geom_ribbon(data = dat.southern_brms.beta.disp,
                            aes(y = response,
                                ymin = lower.HPD, ymax = upper.HPD,
                                x = as.numeric(as.character(Year)),
                                fill = 'brms'),
                            alpha=0.3) +
                geom_line(data = dat.southern_glmmTMB.beta.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'glmmTMB')) +
                geom_ribbon(data = dat.southern_glmmTMB.beta.disp,
                            aes(y = response,
                                ymin = lower.CL, ymax = upper.CL,
                                x = as.numeric(as.character(Year)),
                                fill = 'glmmTMB'),
                            alpha=0.3) + 
                geom_line(data = dat.southern_inla.beta.disp,
                          aes(y = mean,
                              x = as.numeric(as.character(Year)),
                              colour = 'inla')) +
                geom_ribbon(data = dat.southern_inla.beta.disp,
                            aes(y = mean,
                                ymin = lower, ymax = upper,
                                x = as.numeric(as.character(Year)),
                                fill = 'inla'),
                            alpha=0.3) 
            print(g1)
            dev.off()
            }
            ## ----end
            ## ---- beta.ry.disp
            {
            load(file=paste0('../data/modelled/mod.southern_glmmTMB.beta.ry.disp.RData'))
            load(file=paste0('../data/modelled/mod.southern_inla.beta.ry.disp.RData'))
            pdf(file = '../output/figures/comparison.southern.ry.pdf')
            g1 <- ggplot() +
                ## geom_line(data = dat.southern_brms.beta.ry.disp,
                ##           aes(y = response,
                ##               x = as.numeric(as.character(Year)),
                ##               colour = 'brms')) +
                ## geom_ribbon(data = dat.southern_brms.beta.ry.disp,
                ##             aes(y = response,
                ##                 ymin = lower.HPD, ymax = upper.HPD,
                ##                 x = as.numeric(as.character(Year)),
                ##                 fill = 'brms'),
                ##             alpha=0.3) +
                geom_line(data = dat.southern_glmmTMB.beta.ry.disp,
                          aes(y = response,
                              x = as.numeric(as.character(Year)),
                              colour = 'glmmTMB')) +
                geom_ribbon(data = dat.southern_glmmTMB.beta.ry.disp,
                            aes(y = response,
                                ymin = lower.CL, ymax = upper.CL,
                                x = as.numeric(as.character(Year)),
                                fill = 'glmmTMB'),
                            alpha=0.3) + 
                geom_line(data = dat.southern_inla.beta.ry.disp,
                          aes(y = mean,
                              x = as.numeric(as.character(Year)),
                              colour = 'inla')) +
                geom_ribbon(data = dat.southern_inla.beta.ry.disp,
                            aes(y = mean,
                                ymin = lower, ymax = upper,
                                x = as.numeric(as.character(Year)),
                                fill = 'inla'),
                            alpha=0.3) 
            print(g1)
            dev.off()
            }
            ## ----end

        }
    }

    if ( 1 == 2) {
        ## Compare the models old
    {
        ## ---- Southern Compare models


        load(file='../data/modelled/dat.southern.RData')
        dat.southern.original <- dat.southern
        load(file='../data/modelled/mod.southern.RData')
        load(file='../data/modelled/mod.southern_inla_beta.RData')
        load(file='../data/modelled/newdata.southern_inla_beta.RData')
        newdata.southern_beta <- newdata.southern
        load(file='../data/modelled/mod.southern_inla_binomial.RData')
        load(file='../data/modelled/newdata.southern_inla_binomial.RData')
        newdata.southern_binomial <- newdata.southern
        load(file='../data/modelled/dat.southern_inla_tow.RData')
        load(file='../data/modelled/dat.southern_mgcv.RData')
        load(file='../data/modelled/dat.southern_brms.cumulative.RData')
        load(file='../data/modelled/dat.southern_clmm.RData')
        load(file='../data/modelled/dat.southern_brms.beta.RData')
        load(file='../data/modelled/dat.southern_glmmTMB.RData')
        load(file='../data/modelled/dat.southern_mgcv.beta.RData')

        ## original
        g1 <- dat.southern.original %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('Original (stan binomial reef level)') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g1

        ## inla (reef level binomial)
        g2 <- newdata.southern_binomial %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA reef level binomial') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g2

        ## inla (reef level beta)
        g3 <- newdata.southern_beta %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA reef level beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g3

        ## inla (tow level beta)
        g4 <- dat.southern %>%
            ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
            geom_line(color='blue') +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('INLA tow level beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g4

        ## mgcv ordinal
        g5 <- dat.southern.mgcv %>%
            ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('mgcv ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))

        ## brms ordinal
        g6 <- dat.southern_brms.cumulative %>%
            ggplot(aes(y=estimate, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill='blue', alpha=0.3) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g6    

        ## clmm ordinal
        g7 <- dat.southern_clmm %>%
            ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('clmm ordinal') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g7    

        ## brms beta
        g8 <- dat.southern_brms.beta %>%
            ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
            geom_line(color='blue') +  
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_discrete('Raw data aggregate') + 
            rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))
        g8

        library(patchwork)
        g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8


        dat.southern_brms.beta %>%
            ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
            ## geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill='BRMS beta'), alpha=0.3) +
            geom_line(aes(color='BRMS beta')) +
            ## geom_ribbon(data=dat.southern_brms.cumulative, aes(y=estimate, ymin=conf.low, ymax=conf.high, fill='BRMS ordinal'), alpha=0.3) +
            geom_line(data=dat.southern_brms.cumulative, aes(y=estimate, color='BRMS ordinal')) +
            ## geom_ribbon(data=newdata.southern_beta, aes(y=mean, ymin=lower, ymax=upper, fill='INLA beta'), alpha=0.3) +
            geom_line(data=newdata.southern_beta, aes(y=mean, color='INLA beta')) +
            geom_line(data=dat.southern_clmm, aes(y=Mean, color='clmm ordinal')) +
            geom_line(data=dat.southern.mgcv, aes(y=Mean, color='mgcv ordinal')) +
            geom_line(data=dat.southern_glmmTMB, aes(y=response, color='glmmTMB beta')) +
            geom_line(data=dat.southern_mgcv.beta, aes(y=response, color='mgcv beta')) +
            scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
            scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
            scale_color_brewer('Model', type='qual') + 
            scale_fill_discrete('Model') + 
            ## rawAdd +
            theme_classic() +
            theme(axis.title.x=element_blank(),
                  legend.position=c(0.01,0.01), legend.justification=c(0,0),
                  panel.grid.minor=element_line(),
                  panel.grid.major=element_line()) +
            ggtitle('brms beta') +
            guides(color=guide_legend(nrow=2, byrow=TRUE))


        rm(list=c('dat.southern','mod.southern','mod.southern_inla_beta','newdata.southern','newdata.southern_beta', 'mod_southern_inla_binomial','newdata.southern_binomial'))
    }
        ## ----end
        ## ---- junk
    {
        ## summary(mod.southern)
        ## dat.southern = data.frame(Location='Southern',Year=unique(dat.all.southern$Year), N=length(unique(dat.all.southern$REEF_NAME)))
        ## Xmat = model.matrix(~Year, dat.southern)
        ## coefs = data.frame(mod.southern) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
        ## coefs = mod.southern$fit %>% as.data.frame() %>% dplyr::select(matches('b_.*'))
        ## Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
        ## dat.southern = cbind(dat.southern,
        ##             plyr:::adply(Fit,2,function(x) {
        ##     data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
        ##             })
        ## )


        ##     ## glmmTMB
        ##     library(glmmTMB)
        ##     library(emmeans)
        ##     mod.glmmTMB <-  glmmTMB(Cover ~ Year+(Year|REEF_NAME),
        ##                        data=dat.all.southern,
        ##                        #weights=dat.all.gbr$Weight,
        ##                        family=beta_family)
        ##     emmeans(mod.glmmTMB, ~Year, type='response') %>% as.data.frame %>%
        ##     ggplot() +
        ##     geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, x=as.numeric(as.character(Year))), fill='blue', alpha=0.2) +
        ##     geom_line(aes(y=response, x=as.numeric(as.character(Year)))) +
        ##     theme_bw()


        ## dat.southern %>% arrange(Year)
        ## dat.southern %>%
        ##     ggplot(aes(y=mean, x=as.numeric(Year))) +
        ##     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        ##     geom_line()

        ## save(dat.southern, file='../data/modelled/dat.southern.RData')
        ## save(dat.all.southern, file='../data/modelled/dat.all.southern.RData')
        ## save(mod.southern, file='../data/modelled/mod.southern.RData')
        ## rm(list=c('last_year','mod.southern', 'dat.southern', 'coefs', 'Fit'))
        ## gc()
    }
        ## ----end
    }
## ----end
}

## -----END HERE ----------

if (1==2) {
if (model == 'original') {
  mod.central <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|REEF_NAME),
                             data=dat.all.central,
                             family=binomial,
                             iter=5000,
                             warmup=2500,
                             chains=3,cores=3)
    dat.central <- data.frame(Location='Central GBR',Year=unique(dat.all.central$Year), N=length(unique(dat.all.central$REEF_NAME))) 
    Xmat <- model.matrix(~Year, dat.central)

    coefs = data.frame(mod.central) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
    Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
    dat.central = cbind(dat.central,
                    plyr:::adply(Fit,2,function(x) {
                        data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                    })
                    )
    save(dat.central, file='../data/modelled/dat.central.RData')
    save(mod.central, file='../data/modelled/mod.central.RData')
    rm(list=c('last_year','mod.central'))
    gc()

} else if (model=='stan_glmer') {
    mod.central <- Modelstan_glmer(form=Cover ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                  dat=dat.all.central)
    newdata.central <- cellMeansOriginal(mod.central, dat=dat.all.central, location='Northern') 
    newdata.central %>%
        ggplot(aes(y=mean, x=as.numeric(Year))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line()
    save(mod.central, file='../data/modelled/mod.central_stan_glmer.RData')
    save(newdata.central, file='../data/modelled/newdata.central_stan_glmer.RData')
    rm(list=c('last_year','mod.central', 'newdata.central', 'coefs', 'Fit'))
    gc()
} else if (model == "INLA_tow") {
    manta.tow.central = manta.tow %>%
        filter(Region=='Central GBR') %>%
        droplevels
    dat.inla <- dataINLA(dat=manta.tow.central, level='tow')
    dat = dat.inla[['dat.1']]
    dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
    rawAdd <- ggproto_Raw(dat.cellmeans)
    dat.scale = manta.tow.central %>%
        filter(!is.na(REEF_NAME)) %>%
        droplevels %>%
        group_by(REEF_NAME, REPORT_YEAR) %>%
        summarise(Tows=length(TOW_SEQ_NO)) %>%
        ungroup %>% group_by(REEF_NAME) %>%
        summarise(Tows=max(Tows))
    mod.central <- ModelINLA_beta(form=Cover~Year +
                                  f(P_CODE.mod, model='iid') +
                                  f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                  f(REEF_NAME2, YEAR1, model='iid', scale=dat.scale$Tows),
                                  dat=dat,
                              family='beta',
                              weights=NULL)
    
    dat.central <- cellMeansINLA(mod=mod.central, newdata.hcc=dat.inla[['newdata.hcc']],
                             n.2=dat.inla[['n.2']], FUN=plogis)
    dat.raw = dat %>% group_by(REEF_NAME, Year) %>% summarize(Cover=mean(Cover, na.rm=TRUE)) %>% filter(!is.na(Cover)) %>% droplevels
    dat.central %>%
        ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
       rawAdd +
        ggtitle('GBR INLA tow level')

    save(mod.central, file='../data/modelled/mod.central_inla_tow.RData')
    save(dat.central, file='../data/modelled/dat.central_inla_tow.RData')
    rm(list=c('manta.tow.central','dat.central', 'mod.central'))
    gc()
    
} else if (model=='brms') {
    priors <- prior(normal(0, 5), class = "b") +
        prior(normal(0, 5), class = "Intercept") +
        ## prior(gamma(2, 0.1), class = "sd") +
        prior(gamma(1, 0.5), class = "sd") +
        prior(gamma(0.01, 0.01), class = "phi")
    inits = list(list(phi=list(rgamma(1,0.1,0.1))),
                 list(phi=list(rgamma(1,0.1,0.1))),
                 list(phi=list(rgamma(1,0.1,0.1)))
                 )
    mod.central <- Modelbrms(form=Cover|weights(Tows) ~ Year + (1|REEF_NAME),
                              dat=dat.all.central)
    
    mod.central <- Modelbrms(form=Cover|weights(Tows) ~ Year + (Year|REEF_NAME),
                              dat=dat.all.central)
} else if (model=='INLA_binomial') {
    dat.inla <- dataINLA(dat=dat.all.central)
    mod.central <- ModelINLA_binomial(form=Cvr1~Year+
                                  f(P_CODE.mod, model='iid')+
                                  f(REEF_NAME, model='iid') +
                                  f(REEF_NAME1, Year, model='iid'),
                              dat=dat.inla[['dat.1']])
    newdata.central <- cellMeansINLA(mod=mod.central, newdata.hcc=dat.inla[['newdata.hcc']],
                                 n.2=dat.inla[['n.2']])
    newdata.central %>%
        ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
        rawAdd
    save(mod.central, file='../data/modelled/mod.central_inla_binomial.RData')
    save(newdata.central, file='../data/modelled/newdata.central_inla_binomial.RData')
    rm(list=c('mod.central', 'newdata.central'))
    gc()
    
} else if (model=="INLA_beta") {
    dat.inla <- dataINLA(dat=dat.all.central %>% mutate(W=NA))
    dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
    mod.central <- ModelINLA_beta(form=Cover~Year+
                                  f(P_CODE.mod, model='iid')+
                                  f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                  f(REEF_NAME1, Year, model='iid', scale=dat.scale$Tows),
                              dat=dat.inla[['dat.1']])
    newdata.central <- cellMeansINLA(mod=mod.central, newdata.hcc=dat.inla[['newdata.hcc']],
                                 n.2=dat.inla[['n.2']])
    newdata.central %>%
        ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
        rawAdd
    save(mod.central, file='../data/modelled/mod.central_inla_beta.RData')
    save(newdata.central, file='../data/modelled/newdata.central_inla_beta.RData')
    rm(list=c('mod.central', 'newdata.central'))
    gc()
} else if (model=='glmmTMB') {
    library(glmmTMB)
    library(emmeans)
    mod.central <- glmmTMB(Cover ~ Year + (1|P_CODE.mod/REEF_NAME),
                            data=dat.all.central,
                            weights=dat.all.central$Tows,
                            family=beta_family())
    emmeans(mod.central, ~Year, type='response') %>%
        as.data.frame() %>% 
        ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
        geom_line() +
        rawAdd
        
}

## Compare the models

load(file='../data/modelled/dat.central.RData')
dat.central.original <- dat.central
load(file='../data/modelled/mod.central.RData')
load(file='../data/modelled/mod.central_inla_beta.RData')
load(file='../data/modelled/newdata.central_inla_beta.RData')
newdata.central_beta <- newdata.central
load(file='../data/modelled/mod.central_inla_binomial.RData')
load(file='../data/modelled/newdata.central_inla_binomial.RData')
newdata.central_binomial <- newdata.central
load(file='../data/modelled/dat.central_inla_tow.RData')

## original
g1 <- dat.central.original %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('Original (stan binomial reef level)') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g1

## inla (reef level binomial)
g2 <- newdata.central_binomial %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('INLA reef level binomial') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g2

## inla (reef level beta)
g3 <- newdata.central_beta %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('INLA reef level beta') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g3

## inla (tow level beta)
g4 <- dat.central %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('INLA tow level beta') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g4

library(patchwork)
g1 + g2 + g3 + g4

rm(list=c('dat.central','mod.central','mod.central_inla_beta','newdata.central','newdata.central_beta', 'mod_central_inla_binomial','newdata.central_binomial'))

## } else {
##     mod.central <- stan_glmer(Cover ~ Year + (Year|P_CODE.mod/REEF_NAME),
##                           data=dat.all.central,
##                           family=mgcv::betar,
##                           iter=5000,
##                           warmup=2500,
##                           chains=3,cores=3,
##                           adapt_delta=0.95
##                           )
## }
## dat.central = data.frame(Location='Central',Year=unique(dat.all.central$Year), N=length(unique(dat.all.central$REEF_NAME))) 
## Xmat = model.matrix(~Year, dat.central)
## coefs = data.frame(mod.central) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
## Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
## dat.central = cbind(dat.central,
##             plyr:::adply(Fit,2,function(x) {
##     data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
## })
## )
## save(dat.central, file='../data/modelled/dat.central.RData')
## save(dat.all.central, file='../data/modelled/dat.all.central.RData')
## save(mod.central, file='../data/modelled/mod.central.RData')
## rm(list=c('last_year','mod.central', 'dat.central', 'coefs', 'Fit'))
## gc()

## Southern

dat.all.southern = dat.all %>%
    filter(Location=='Southern GBR') %>%
    droplevels %>% 
    mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
    group_by(REEF_NAME) %>%
    mutate(W=mean(Tows, na.rm=TRUE)) %>%
    ungroup %>%
    mutate(W1=W/sum(W)) %>%
    group_by(Year) %>%
    mutate(W2=Tows/sum(Tows)) %>%
    ungroup
save(dat.all.southern, file='../data/modelled/dat.all.southern.RData')
## dat.all.southern = dat.all %>% filter(Location=='Southern GBR') %>% droplevels
#dat.all.southern$Cvr1 = as.integer(as.vector(dat.all.southern$Cover) * dat.all.southern$Tows)
#dat.all.southern$Cvr0 = dat.all.southern$Tows - dat.all.southern$Cvr1
dat.all.southern.cellmeans <- cellMeansRaw(dat.all.southern)
rawAdd <- ggproto_Raw(dat.all.southern.cellmeans)

manta.tow.southern.cellmeans <- cellMeansRaw(manta.tow %>%
                                             filter(Region=='Northern GBR') %>%
                                             droplevels %>% 
                                             group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                             summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                             ungroup)
if (model == 'original') {
  mod.southern <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|REEF_NAME),
                              data=dat.all.southern,
                              family=binomial,
                              iter=5000,
                              warmup=2500,
                              chains=3,cores=3)
    dat.southern <- data.frame(Location='Southern GBR',Year=unique(dat.all.southern$Year), N=length(unique(dat.all.southern$REEF_NAME))) 
    Xmat <- model.matrix(~Year, dat.southern)

    coefs = data.frame(mod.southern) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
    Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
    dat.southern = cbind(dat.southern,
                    plyr:::adply(Fit,2,function(x) {
                        data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                    })
                    )
    save(dat.southern, file='../data/modelled/dat.southern.RData')
    save(mod.southern, file='../data/modelled/mod.southern.RData')
    rm(list=c('last_year','mod.southern'))
    gc()

} else if (model=='stan_glmer') {
    mod.southern <- Modelstan_glmer(form=Cover ~ Year+(Year|P_CODE.mod/REEF_NAME),
                                  dat=dat.all.southern)
    newdata.southern <- cellMeansOriginal(mod.southern, dat=dat.all.southern, location='Northern') 
    newdata.southern %>%
        ggplot(aes(y=mean, x=as.numeric(Year))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line()
    save(mod.southern, file='../data/modelled/mod.southern_stan_glmer.RData')
    save(newdata.southern, file='../data/modelled/newdata.southern_stan_glmer.RData')
    rm(list=c('last_year','mod.southern', 'newdata.southern', 'coefs', 'Fit'))
    gc()
} else if (model == "INLA_tow") {
    manta.tow.southern = manta.tow %>%
        filter(Region=='Southern GBR') %>%
        droplevels
    dat.inla <- dataINLA(dat=manta.tow.southern, level='tow')
    dat = dat.inla[['dat.1']]
    dat.cellmeans <- cellMeansRaw(dat %>% mutate(Tows=1) %>% filter(!is.na(Cover)))
    rawAdd <- ggproto_Raw(dat.cellmeans)
    dat.scale = manta.tow.southern %>%
        filter(!is.na(REEF_NAME)) %>%
        droplevels %>%
        group_by(REEF_NAME, REPORT_YEAR) %>%
        summarise(Tows=length(TOW_SEQ_NO)) %>%
        ungroup %>% group_by(REEF_NAME) %>%
        summarise(Tows=max(Tows))
    mod.southern <- ModelINLA_beta(form=Cover~Year +
                                  f(P_CODE.mod, model='iid') +
                                  f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                  f(REEF_NAME2, YEAR1, model='iid', scale=dat.scale$Tows),
                                  dat=dat,
                              family='beta',
                              weights=NULL)
    
    dat.southern <- cellMeansINLA(mod=mod.southern, newdata.hcc=dat.inla[['newdata.hcc']],
                             n.2=dat.inla[['n.2']], FUN=plogis)
    dat.raw = dat %>% group_by(REEF_NAME, Year) %>% summarize(Cover=mean(Cover, na.rm=TRUE)) %>% filter(!is.na(Cover)) %>% droplevels
    dat.southern %>%
        ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
       rawAdd +
        ggtitle('GBR INLA tow level')

    save(mod.southern, file='../data/modelled/mod.southern_inla_tow.RData')
    save(dat.southern, file='../data/modelled/dat.southern_inla_tow.RData')
    rm(list=c('manta.tow.southern','dat.southern', 'mod.southern'))
    gc()
    
} else if (model=='brms') {
    priors <- prior(normal(0, 5), class = "b") +
        prior(normal(0, 5), class = "Intercept") +
        ## prior(gamma(2, 0.1), class = "sd") +
        prior(gamma(1, 0.5), class = "sd") +
        prior(gamma(0.01, 0.01), class = "phi")
    inits = list(list(phi=list(rgamma(1,0.1,0.1))),
                 list(phi=list(rgamma(1,0.1,0.1))),
                 list(phi=list(rgamma(1,0.1,0.1)))
                 )
    mod.southern <- Modelbrms(form=Cover|weights(Tows) ~ Year + (1|REEF_NAME),
                              dat=dat.all.southern)
    
    mod.southern <- Modelbrms(form=Cover|weights(Tows) ~ Year + (Year|REEF_NAME),
                              dat=dat.all.southern)
} else if (model=='INLA_binomial') {
    dat.inla <- dataINLA(dat=dat.all.southern)
    mod.southern <- ModelINLA_binomial(form=Cvr1~Year+
                                  f(P_CODE.mod, model='iid')+
                                  f(REEF_NAME, model='iid') +
                                  f(REEF_NAME1, Year, model='iid'),
                              dat=dat.inla[['dat.1']])
    newdata.southern <- cellMeansINLA(mod=mod.southern, newdata.hcc=dat.inla[['newdata.hcc']],
                                 n.2=dat.inla[['n.2']])
    newdata.southern %>%
        ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
        rawAdd
    save(mod.southern, file='../data/modelled/mod.southern_inla_binomial.RData')
    save(newdata.southern, file='../data/modelled/newdata.southern_inla_binomial.RData')
    rm(list=c('mod.southern', 'newdata.southern'))
    gc()
    
} else if (model=="INLA_beta") {
    dat.inla <- dataINLA(dat=dat.all.southern %>% mutate(W=NA))
    dat.scale = dat.inla[['dat.1']] %>% filter(!is.na(REEF_NAME)) %>% droplevels %>% group_by(REEF_NAME) %>% summarise(Tows=max(Tows)) 
    mod.southern <- ModelINLA_beta(form=Cover~Year+
                                  f(P_CODE.mod, model='iid')+
                                  f(REEF_NAME, model='iid', scale=dat.scale$Tows) +
                                  f(REEF_NAME1, Year, model='iid', scale=dat.scale$Tows),
                              dat=dat.inla[['dat.1']])
    newdata.southern <- cellMeansINLA(mod=mod.southern, newdata.hcc=dat.inla[['newdata.hcc']],
                                 n.2=dat.inla[['n.2']])
    newdata.southern %>%
        ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
        rawAdd
    save(mod.southern, file='../data/modelled/mod.southern_inla_beta.RData')
    save(newdata.southern, file='../data/modelled/newdata.southern_inla_beta.RData')
    rm(list=c('mod.southern', 'newdata.southern'))
    gc()
} else if (model=='glmmTMB') {
    library(glmmTMB)
    library(emmeans)
    mod.southern <- glmmTMB(Cover ~ Year + (1|P_CODE.mod/REEF_NAME),
                            data=dat.all.southern,
                            weights=dat.all.southern$Tows,
                            family=beta_family())
    emmeans(mod.southern, ~Year, type='response') %>%
        as.data.frame() %>% 
        ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
        geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
        geom_line() +
        rawAdd
        
}


## Compare the models

load(file='../data/modelled/dat.southern.RData')
dat.southern.original <- dat.southern
load(file='../data/modelled/mod.southern.RData')
load(file='../data/modelled/mod.southern_inla_beta.RData')
load(file='../data/modelled/newdata.southern_inla_beta.RData')
newdata.southern_beta <- newdata.southern
load(file='../data/modelled/mod.southern_inla_binomial.RData')
load(file='../data/modelled/newdata.southern_inla_binomial.RData')
newdata.southern_binomial <- newdata.southern
load(file='../data/modelled/dat.southern_inla_tow.RData')

## original
g1 <- dat.southern.original %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('Original (stan binomial reef level)') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g1

## inla (reef level binomial)
g2 <- newdata.southern_binomial %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('INLA reef level binomial') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g2

## inla (reef level beta)
g3 <- newdata.southern_beta %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('INLA reef level beta') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g3

## inla (tow level beta)
g4 <- dat.southern %>%
    ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) + 
    scale_color_discrete('Raw data aggregate') + 
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('INLA tow level beta') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
g4

library(patchwork)
g1 + g2 + g3 + g4

rm(list=c('dat.southern','mod.southern','mod.southern_inla_beta','newdata.southern','newdata.southern_beta', 'mod_southern_inla_binomial','newdata.southern_binomial'))


## } else if (model=='stan_glmer') {
##     mod.southern <- stan_glmer(Cover ~ Year + (Year|P_CODE.mod/REEF_NAME),
##                           data=dat.all.southern,
##                           family=mgcv::betar,
##                           iter=5000,
##                           warmup=2500,
##                           chains=3,cores=3,
##                           adapt_delta=0.95
##                           )
## }

## dat.southern = data.frame(Location='Southern',Year=unique(dat.all.southern$Year), N=length(unique(dat.all.southern$REEF_NAME))) 
## Xmat = model.matrix(~Year, dat.southern)
## coefs = data.frame(mod.southern) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
## Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
## dat.southern = cbind(dat.southern,
##             plyr:::adply(Fit,2,function(x) {
##     data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
## })
## )
## save(dat.southern, file='../data/modelled/dat.southern.RData')
## save(dat.all.southern, file='../data/modelled/dat.all.southern.RData')
## save(mod.southern, file='../data/modelled/mod.southern.RData')
## rm(list=c('last_year','mod.southern', 'dat.southern', 'coefs', 'Fit'))
## gc()





## There may be a need to produce individual reef trend estimates ===========================================
## Northern reefs


MONITORING_model_spatiotemporal_fig <- function(newdata, region) {
    ## This is modified from ALT project
    load('data/spatial/gbr_3Zone.RData')
    gbr_3Zone <- sf::st_as_sf(whagbr.n + whagbr.c + whagbr.s) %>%
        sf::st_set_crs('EPSG:4283') %>%
        mutate(Region=c('Northern GBR','Central GBR','Southern GBR'))
    load('data/spatial/qld.RData')
    gbr <- gbr_3Zone %>% filter(if(region!='GBR') Region==region else Region!='GBR')
    bb <- sf::st_bbox(gbr) + c(-1, -1, 1, 1)
    g <- newdata %>% ungroup %>%
        sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(gbr)) %>%
        ggplot() +
        sf::geom_sf(data = qld %>% st_crop(bb), fill = "grey", color = "grey40") +
        sf::geom_sf(aes(color = Mean)) +
        scale_color_gradientn('Coral Cover', colours=grDevices::blues9[-1])+
        #annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
        #annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
        #annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
        #annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
        #annotation_scale(location = "tr", pad_y=unit(2, units='cm')) +
        sf::geom_sf(data = gbr, fill = NA) +
        sf::coord_sf(expand=FALSE) +
        theme_bw() +
        theme(#legend.position = c(1,0),legend.justification=c(1,0),
            legend.background = element_rect(color='black', fill='white'),
            axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text = element_text(size=6)) +
        facet_wrap(~Year, ncol=6)
    return(g)
}


MONITORING_weighted_model_temporal_fig <- function(fit, dat, full.dat, region) {
    d = MONITORING_weighted_compile(fit=fit, dat=dat, full.dat=full.dat, region=region) %>%
        dplyr::select(-Cell) %>%
        mutate(Value=Value*Weights) %>%
        group_by(Sample, Year) %>%
        summarize(Value=binomial()$linkinv(sum(Value)),
                  N=mean(N)) %>%
        ungroup %>%
        group_by(Year)
    d1 = d %>%
        tidybayes::point_intervalh(Value, .point=mean)
    d = d %>% summarize(N=mean(N)) %>%
        left_join(d1) %>%
        dplyr::rename(mean=Value, lower=.lower, upper=.upper)
    g = MONITORING_model_fig(d)
    return(list(d=d, g=g))
}


MONITORING_weighted_compile <- function(fit, dat, full.dat, region) {
    ## generate new prediction dataset (year/reef level)
    nd = expand.grid(Year=factor(levels(dat$Year)),
                     REEF_NAME=factor(sort(unique(dat$REEF_NAME))),
                     N=length(unique(dat$REEF_NAME))) %>%
        left_join(dat %>% ungroup %>%
                  group_by(REEF_NAME) %>%
                  summarize(Tows=mean(Tows, na.rm=TRUE)) %>%
                  ungroup) %>%
        mutate(REEF_NAME = factor(REEF_NAME, levels=unique(sort(dat$REEF_NAME))))
    ## generate reef weights (based on average number of tows over time)
    wts = dat %>% dplyr::select(REEF_NAME, Tows) %>%
        group_by(REEF_NAME) %>%
        summarize(Tows=mean(Tows, na.rm=TRUE)) %>%
        ungroup %>%
        mutate(Weights = Tows/sum(Tows))
  nms <- rownames(ranef(fit)[[1]])
  cols <- colnames(ranef(fit)[[1]])[-1]
  nd = expand.grid(Year=factor(cols), REEF_NAME=factor(nms))

    Xmat = model.matrix(~Year, data=nd) ## May need to remove Latitude*Longitude
    Xmat1 = model.matrix(~0+Year:REEF_NAME, data=nd)
    Xmat = cbind(Xmat, Xmat1)
  ## Extract the relevant coefficients
  ## reg <- paste0('^X.Intercept.|^Year.*|^b.*',nms,'.*', collapse='|')
  reg <- paste0('^X.Intercept.|^Year.*|^b.*')
    ## coefs = data.frame(fit) %>% dplyr:::select(matches('^X.Intercept.|^Year.*|^b.*'))
  coefs = data.frame(fit) %>% dplyr:::select(matches(reg))
  coefs = coefs[1:ncol(Xmat), 1:ncol(Xmat)]
    ## Perform matrix multiplications
    Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))  #May need to move the back transform to after weighting and allow other link functions
    ## Add the the MCMC summaries back to the newdata (not this is still on link scale
    ## so that we can perform weighted aggregations)
    nd = cbind(nd,
                    plyr:::adply(Fit,2,function(x) {
                        data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                    })
               )
    reef = "ASHMORE BANKS (3)"
    ggplot(nd %>% filter(REEF_NAME==reef)) +
        geom_point(aes(y=mean, x=Year)) + geom_ribbon(aes(ymin=lower, ymax=upper, x=as.numeric(Year)), fill='blue', alpha=0.3) + geom_line(aes(y=mean, x=as.numeric(Year))) +
        geom_point(data=dat %>% filter(REEF_NAME==reef), aes(y=Cover, x=Year), color='red') +
        geom_line(data=dat %>% filter(REEF_NAME==reef), aes(y=Cover, x=as.numeric(Year)), color='red')
    
    ##Return the full MCMC sample matrix joined to the newdata and weights
    ##NOTE, this is all on the link scale so that weighting can occur on the link scale
    return(Fit %>% as.data.frame %>%
           mutate(Sample=1:n()) %>% 
           gather(key=Cell,value=Value, -Sample) %>%
           mutate(Cell=as.numeric(as.character(Cell))) %>% 
           full_join(nd %>% mutate(Cell = 1:n())) %>%
           left_join(wts %>% dplyr::select(REEF_NAME, Weights))
           )
}


load(file='../data/modelled/mod.northern.RData')
load(file='../data/modelled/dat.all.northern.RData')
cellmeans = MONITORING_weighted_model_temporal_fig(fit = mod.northern,

                                                   full.dat=manta.sum,
                                                   region='Northern GBR')
}
