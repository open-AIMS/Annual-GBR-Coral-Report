source('CoralTrends_functions.R') 
CoralTrends_checkPackages()
source('CoralTrends_config.R')

## ----PrepareData
load('data/processed/manta.sum.RData')
library(xtable)

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
    dplyr:::select(Cover, REEF_NAME, Tows,P_CODE.mod,Location,REPORT_YEAR) %>%
    mutate(Year=factor(REPORT_YEAR), N=length(unique(REEF_NAME))) %>% ungroup() %>%
    mutate(Cvr1 = as.integer(as.vector(Cover) * Tows), Cvr0 = Tows - Cvr1)



## Calculate some relatively simple comparisons
dat.all.gbr = dat.all %>% droplevels
## ----


## ----DefinePrePost
dat.all.gbr = dat.all.gbr %>%
    mutate(Comp='None',
           Comp = ifelse(Location=='Northern' & REPORT_YEAR==coralchange_northern.y1, 'Pre', Comp),
           Comp = ifelse(Location=='Northern' & REPORT_YEAR==coralchange_northern.y2, 'Post', Comp),
           Comp = ifelse(Location=='Central' & REPORT_YEAR==coralchange_central.y1, 'Pre', Comp),
           Comp = ifelse(Location=='Central' & REPORT_YEAR==coralchange_central.y2, 'Post', Comp),
           Comp = ifelse(Location=='Southern' & REPORT_YEAR==coralchange_southern.y1, 'Pre', Comp),
           Comp = ifelse(Location=='Southern' & REPORT_YEAR==coralchange_southern.y2, 'Post', Comp)
           ) %>% filter(Comp!='None') %>%
    mutate(Location=factor(Location, levels=c('Northern','Central','Southern')))
## ----
## If we aggregate over reefs first

## ----PrePostEDA
dat.all.gbr %>% ggplot(aes(x=Cover, fill=Comp)) +
    geom_density(alpha=0.3) +
    facet_grid(~Location) +
    scale_fill_discrete('') +
    theme_bw()
## ----

## ----PrePostApproach1
dat.all.gbr.sum=(dat.all.gbr.cast = dat.all.gbr %>% group_by(Location, Comp) %>%
    summarize(Cover=mean(Cover, na.rm=TRUE), Tows=sum(Tows), N=n()) %>%
    ungroup %>%
        group_by(Location) %>%
        mutate(Tows=mean(Tows), N=mean(N)) %>%
        spread(key=Comp, value=Cover)) %>%
        ungroup %>%
    mutate(Diff=Pre-Post, DiffP=100*Diff/Pre)
dat.all.gbr.sum = dat.all.gbr.sum %>%
    full_join(dat.all.gbr.sum %>% summarize(DiffP=weighted.mean(DiffP,w=Tows)) %>% mutate(Location='GBR'))
print(xtable(dat.all.gbr.sum),booktabs=TRUE,comment=FALSE,include.rownames=FALSE, hline=c(-1,0,3,nrow(dat.all.gbr.sum)))
## ----
## ----PrePostApproach1.sum
dat.all.gbr.sum %>% summarize(DiffP=weighted.mean(DiffP,w=Tows))
## ----

## Median based
## ----PrePostApproach2
dat.all.gbr.sum=(dat.all.gbr.cast = dat.all.gbr %>% group_by(Location, Comp) %>%
    summarize(Cover=median(Cover, na.rm=TRUE), Tows=sum(Tows), N=n()) %>%
    ungroup %>%
        group_by(Location) %>%
        mutate(Tows=median(Tows), N=mean(N)) %>%
        spread(key=Comp, value=Cover)) %>%
        ungroup %>%
    mutate(Diff=Pre-Post, DiffP=100*Diff/Pre)
dat.all.gbr.sum = dat.all.gbr.sum %>%
    full_join(dat.all.gbr.sum %>% summarize(DiffP=spatstat::weighted.median(DiffP,w=Tows)) %>% mutate(Location='GBR'))
print(xtable(dat.all.gbr.sum),booktabs=TRUE,comment=FALSE,include.rownames=FALSE, hline=c(-1,0,3,nrow(dat.all.gbr.sum)))
## ----
dat.all.gbr.sum
dat.all.gbr.sum %>% summarize(DiffP=spatstat::weighted.median(DiffP,w=Tows))



## If we calculate the difference per reef first
## ---- PrePostApproach3
dat.all.gbr.cast.1=dat.all.gbr %>% group_by(Location, REEF_NAME) %>%
     mutate(Tows=mean(Tows,na.rm=TRUE)) %>%
     dplyr::select(Location,REEF_NAME,Tows,Comp,Cover) %>%
     spread(key=Comp, value=Cover) %>%
     ungroup %>%
     mutate(Diff=Pre-Post, DiffP=100*Diff/Pre)
#dat.all.gbr.cast.1
dat.all.gbr.cast.1 %>% ggplot(aes(x=DiffP)) + geom_density(fill='blue', alpha=0.3) + facet_grid(~Location) + theme_bw()
dat.all.gbr.sum.1=dat.all.gbr.cast.1 %>% na.omit %>% group_by(Location) %>%
    summarize(DiffP=weighted.mean(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows), N=n())
dat.all.gbr.sum.2=dat.all.gbr.sum.1 %>% na.omit%>%
    summarize(DiffP=weighted.mean(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows), N=n())%>% mutate(Location='GBR')
dat.all.gbr.sum.3=dat.all.gbr.cast.1 %>% na.omit %>%
    summarize(DiffP=weighted.mean(DiffP, w=Tows, na.rm=TRUE),Tows=sum(Tows),N=n()) %>% mutate(Location='GBR total')
dat.all.gbr.sum = dat.all.gbr.sum.1 %>% full_join(dat.all.gbr.sum.2) %>% full_join(dat.all.gbr.sum.3)
print(xtable(dat.all.gbr.sum %>% dplyr::select(Location,Tows,N,DiffP)),
      booktabs=TRUE, comment=FALSE, include.rownames=FALSE, hline=c(-1,0,3, 4,nrow(dat.all.gbr.sum)))
## ----
dat.all.gbr.sum.1=dat.all.gbr.cast.1 %>% filter(!is.na(Pre), !is.na(Post)) %>% 
    group_by(Location) %>% summarize(DiffP=weighted.mean(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows))
dat.all.gbr.sum.2=dat.all.gbr.sum.1 %>% summarize(DiffP=weighted.mean(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows))%>% mutate(Location='GBR')
dat.all.gbr.sum.2


## ---- PrePostApproach4
dat.all.gbr.cast.1=dat.all.gbr %>% group_by(Location, REEF_NAME) %>%
     mutate(Tows=median(Tows,na.rm=TRUE)) %>%
     dplyr::select(Location,REEF_NAME,Tows,Comp,Cover) %>%
     spread(key=Comp, value=Cover) %>%
     ungroup %>%
     mutate(Diff=Pre-Post, DiffP=100*Diff/Pre)
#dat.all.gbr.cast.1
dat.all.gbr.cast.1 %>% ggplot(aes(x=DiffP)) + geom_density(fill='blue', alpha=0.3) + facet_grid(~Location) + theme_bw()
dat.all.gbr.sum.1=dat.all.gbr.cast.1 %>% na.omit %>% group_by(Location) %>%
    summarize(DiffP=spatstat::weighted.median(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows))
dat.all.gbr.sum.2=dat.all.gbr.sum.1 %>% na.omit%>%
    summarize(DiffP=spatstat::weighted.median(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows))%>% mutate(Location='GBR')
dat.all.gbr.sum.3=dat.all.gbr.cast.1 %>% na.omit %>%
    summarize(DiffP=spatstat::weighted.median(DiffP, w=Tows, na.rm=TRUE),Tows=sum(Tows)) %>% mutate(Location='GBR total')
dat.all.gbr.sum = dat.all.gbr.sum.1 %>% full_join(dat.all.gbr.sum.2) %>% full_join(dat.all.gbr.sum.3)
print(xtable(dat.all.gbr.sum %>% dplyr::select(Location,Tows,DiffP)),
      booktabs=TRUE, comment=FALSE, include.rownames=FALSE, hline=c(-1,0,3, 4,nrow(dat.all.gbr.sum)))
## ----
dat.all.gbr.sum.1=dat.all.gbr.cast.1 %>% filter(!is.na(Pre), !is.na(Post)) %>% 
    group_by(Location) %>% summarize(DiffP=spatstat::weighted.median(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows))
dat.all.gbr.sum.2=dat.all.gbr.sum.1 %>% summarize(DiffP=spatstat::weighted.median(DiffP, w=Tows, na.rm=TRUE), Tows=sum(Tows))%>% mutate(Location='GBR')
dat.all.gbr.sum.2


## ---- PrePostApproach5
library(nlme)
library(lsmeans)
a = lme(Cover ~ Comp*Location, random=~1|REEF_NAME, data=dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Pre','Post'))))
cc1=print(lsmeans(a, pairwise ~ Comp|Location))
fit=cc1[[1]][1:2,3]
fit.1=(fit[1]-fit[2])/fit[1]
fit=cc1[[1]][3:4,3]
fit.2=(fit[1]-fit[2])/fit[1]
fit=cc1[[1]][5:6,3]
fit.3=(fit[1]-fit[2])/fit[1]
cc2=print(lsmeans(a, pairwise ~ Comp))
fit=cc2[[1]][1:2,2]
fit.4=(fit[1]-fit[2])/fit[1]
## ----
## ----PrePostApproach5a
dat=data.frame(Location=c('Northern','Central','Southern','GBR'),
               DiffP=c(fit.1, fit.2,fit.3,fit.4))
print(xtable(dat),
      booktabs=TRUE, comment=FALSE, include.rownames=FALSE, hline=c(-1,0,3,nrow(dat)))
## ----

library(MASS)
aa = glmmPQL(Cover ~ Comp*Location, random=~1|REEF_NAME, data=dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Pre','Post'))), family=binomial(), weights=Tows)
cc=print(lsmeans(aa, pairwise ~ Comp|Location), type='response')
fit=cc[[1]][1:2,3]
(fit[1]-fit[2])/fit[1]
fit=cc[[1]][3:4,3]
(fit[1]-fit[2])/fit[1]
fit=cc[[1]][5:6,3]
(fit[1]-fit[2])/fit[1]

cc=print(lsmeans(aa, pairwise ~ Comp), type='response')
fit=cc[[1]][1:2,2]
(fit[1]-fit[2])/fit[1]



## lsmeans(aa, 'Comp', by='Location', type='response')

## lsmeans(aa, 'Comp', by='Location', type='response') %>% regrid %>% pairs
## lsmeans(aa, 'Comp', by='Location') %>% regrid %>% pairs(type='response')
## lsmeans(aa, 'Comp', by='Location') %>% regrid(transform='log') %>% pairs()
## lsmeans(aa, 'Comp', by='Location') %>% regrid(transform='none') %>% pairs()
## lsmeans(aa, specs=pairwise~Comp|Location,type='response')

## lsmeans(aa, specs=pairwise~Comp,type='response')

## lsmeans(aa, specs=~1|Location,type='response', weights='equal')
## lsmeans(aa, specs=~1|Location,type='link', weights='equal')
## lsmeans(aa, specs=~Comp,type='response', weights=c(1,2,3))
## lsmeans(aa, specs=pairwise~Comp,type='response')
## lsmeans(aa, specs=pairwise~Comp)

## ---- PrePostApproach6
library(lme4)
bb = glmer(Cover ~ Comp*Location + (1|REEF_NAME), data=dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Pre','Post'))),
           family=binomial(),
           weights=Tows)

cc=print(lsmeans(bb, specs=pairwise~Comp|Location), type='response')
fit=cc[[1]][1:2,3]
fit.1=(fit[1]-fit[2])/fit[1]

fit=cc[[1]][3:4,3]
fit.2=(fit[1]-fit[2])/fit[1]

fit=cc[[1]][5:6,3]
fit.3=(fit[1]-fit[2])/fit[1]

cc=print(lsmeans(bb, specs=pairwise~Comp), type='response')
fit=cc[[1]][1:2,2]
fit.4=(fit[1]-fit[2])/fit[1]
## ----
## ----PrePostApproach6a
dat=data.frame(Location=c('Northern','Central','Southern','GBR'),
               DiffP=c(fit.1, fit.2,fit.3,fit.4))
print(xtable(dat),
      booktabs=TRUE, comment=FALSE, include.rownames=FALSE, hline=c(-1,0,3,nrow(dat)))
## ----

## lsmeans(bb, 'Comp', by='Location', type='response')
## lsmeans(bb, specs=pairwise~Comp,type='response')
## lsmeans(bb, 'Comp', by='Location') %>% regrid %>% pairs(type='response')
## lsmeans(bb, 'Comp') %>% regrid %>% pairs(type='response')

## ---- PrePostApproach7
library(glmmTMB)
bb = glmmTMB(Cover ~ Comp*Location + (1|REEF_NAME), family=list(family='beta', link='logit'),
             data=dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Post','Pre'))), weights=Tows)

recover.data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    recover.data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}

lsm.basis.glmmTMB <- function (object, trms, xlev, grid, vcov.,
                               mode = "asymptotic", component="cond", ...) {
    if (mode != "asymptotic") stop("only asymptotic mode is available")
    if (component != "cond") warning("only tested for conditional component")
    if (missing(vcov.)) 
        V <- as.matrix(vcov(object)[[component]])
    else V <- as.matrix(.my.vcov(object, vcov.))
    dfargs = misc = list()
    if (mode == "asymptotic") {
        dffun = function(k, dfargs) NA
    }
    ## use this? misc = .std.link.labels(family(object), misc)
    ## (used to populate the reminder 
    contrasts = attr(model.matrix(object), "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = fixef(object)[[component]]
    if (length(bhat) < ncol(X)) {
        kept = match(names(bhat), dimnames(X)[[2]])
        bhat = NA * X[1, ]
        bhat[kept] = fixef(object)[[component]]
        modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
        nbasis = estimability::nonest.basis(modmat)
    }
    else nbasis = estimability::all.estble
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
        dfargs = dfargs, misc = misc)
}

cc=print(lsmeans(bb, specs=pairwise~Comp|Location))
fit=binomial()$linkinv(cc[[1]][1:2,3])
fit.1=(fit[2]-fit[1])/fit[2]

fit=binomial()$linkinv(cc[[1]][3:4,3])
fit.2=(fit[2]-fit[1])/fit[2]

fit=binomial()$linkinv(cc[[1]][5:6,3])
fit.3=(fit[2]-fit[1])/fit[2]

cc=print(lsmeans(bb, specs=pairwise~Comp))
fit=binomial()$linkinv(cc[[1]][1:2,2])
fit.4=(fit[2]-fit[1])/fit[2]
## ----
## ----PrePostApproach7a
dat=data.frame(Location=c('Northern','Central','Southern','GBR'),
               DiffP=c(fit.1, fit.2,fit.3,fit.4))
print(xtable(dat),
      booktabs=TRUE, comment=FALSE, include.rownames=FALSE, hline=c(-1,0,3,nrow(dat)))
## ----

newdata = dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Post','Pre')))
Xmat = model.matrix(~Comp*Location, data=newdata)

Xmat = Xmat %>% cbind(newdata %>% dplyr::select(Comp,Location,Tows)) %>%
    group_by(Comp) %>% summarize_if(is.numeric,funs(weighted.mean(., w=Tows,na.rm=TRUE))) %>% dplyr::select(-Comp,-Tows) %>%
    as.matrix

#Xmat = Xmat %>% cbind(newdata %>% dplyr::select(Comp,Location,Tows)) %>%
#    group_by(Comp) %>% summarize_if(is.numeric,funs(mean(., na.rm=TRUE))) %>% dplyr::select(-Comp,-Tows) %>%
#    as.matrix

coefs=fixef(bb)[[1]]
(fit=binomial()$linkinv(coefs %*% t(Xmat)))
(fit[2]-fit[1])/fit[2]


newdata = with(dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Post','Pre'))), expand.grid(Comp=levels(Comp), Location=levels(Location)))
Xmat = model.matrix(~Comp*Location, data=newdata)

Xmat = Xmat %>% cbind(newdata %>% dplyr::select(Comp,Location)) %>%
    group_by(Comp) %>% summarize_if(is.numeric,funs(mean(., na.rm=TRUE))) %>% dplyr::select(-Comp) %>%
    as.matrix

coefs=fixef(bb)[[1]]
(fit=binomial()$linkinv(coefs %*% t(Xmat)))
(fit[2]-fit[1])/fit[2]



lsmeans(bb, 'Comp', by='Location', type='response')
lsmeans(bb, 'Comp', by='Location')
fit=binomial()$linkinv(c(-1.8384800,-0.8675317))
(fit[2]-fit[1])/fit[2]
lsmeans(bb, specs=pairwise~Comp,type='response')
lsmeans(bb, specs=pairwise~Comp)
lsmeans(bb, specs=pairwise~Comp|Location)
lsmeans(bb, 'Comp', by='Location') %>% regrid %>% pairs(type='response')
lsmeans(bb, 'Comp') %>% regrid %>% pairs(type='response')


#library(brms)
#bb = brm(Cover ~ Comp*Location + (1|REEF_NAME), data=dat.all.gbr %>% mutate(Comp=factor(Comp, levels=c('Post','Pre'))),
#           family=binomial(),
#           weights=dat.all.gbr$Tows)
## ----
