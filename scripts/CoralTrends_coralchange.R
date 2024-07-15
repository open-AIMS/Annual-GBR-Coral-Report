library(tidyverse)
library(INLA)

missing = LTMP_checkPackages()

DATA_PATH <<- "../data/"
RDATA_FILE <<- "reef_data.RData"
DATA_FROM='FILE'
FILENAME <<- "data-manta-2024.csv" #"reef_data.csv" #'GBR/LTMP-PT-COVER-GBR.csv'
DOMAIN_NAME <<- "reef" #'GBR'
DATA_TYPE <<- "COVER"
DATA_PROGRAM <<- "LTMP"
CSV_FILE <<- FILENAME
ZIP_FILE <<- gsub('\\.csv','\\.zip', FILENAME)
RDATA_FILE <<- gsub('\\.csv','\\.RData',FILENAME)
FILENAME <<- gsub('\\.csv', '', FILENAME)
DOMAIN_CATEGORY <<- "reef" 
DOMAIN_NAME <<- gsub('.*process/[a-zA-Z]*/[0-9]{4}/[a-zA-Z]*/[^/]*/([^/]*)/.*','\\1',AWS_PATH)
DATA_PROGRAM <<- gsub('.*process/([^/]*).*','\\1',AWS_PATH)
DATA_METHOD <<- gsub('[^/]*//[^/]*/[^/]*/([^/]*)/.*','\\1',AWS_PATH)
LOG_FILE <<- paste0(DATA_PATH, 'log/', FILENAME, '.log')

## system("scp ssm-user@i-03e2a5b592e7eb26a:~/LTMP_web_reporting/data/processed/reef_data.RData ../data/processed/reef_data_PC_.RData")

## data <- read_csv(paste0(DATA_PATH, "primary/", CSV_FILE), col_types = 'cdccccdddTdcdc', trim_ws=TRUE)
## data <- load(paste0(DATA_PATH, "modelled/ann_report_tow.RData"))
## save(data,  file=paste0(DATA_PATH, 'primary/', RDATA_FILE))



data <- get(load("../data/processed/manta.tow.RData"))

data <- data |> dplyr::mutate(REEF = REEF_NAME)

data.shelf <- data %>% dplyr::select(REEF_NAME, SHELF) %>% distinct()
## data.spatial =
##   data %>%
##   group_by(REEF) %>%
##   summarise(Longitude=mean(LONGITUDE,na.rm=TRUE),
##     Latitude=mean(LATITUDE, na.rm=TRUE)) %>%
##   ungroup %>%
##   LTMP_assignSpatialDomain_Zones() %>%
##   LTMP_assignSpatialDomain_NRMRegions() %>%
##   LTMP_assignSpatialDomain_Bioregions() %>%
##   LTMP_assignSpatialDomain_Sectors() %>%
##   LTMP_assignSpatialDomain_Shelf() %>%
##   mutate(Shelf=ifelse(Shelf %in% c('Enclosed Coastal','Open Coastal','Midshelf'), 'Inshore', 'Offshore')) %>%
##   LTMP_assignSpatialDomain_Shelf_from_database(data.shelf) %>%  #this overwrites the previous
##   suppressMessages() %>%
##   suppressWarnings()
## save(data.spatial,  file=paste0(DATA_PATH, 'processed/', gsub('\\.RData', '.spatial.RData', RDATA_FILE)))



LTMP_calcPercent = function(x) {
    ifelse(x=='0', 0,
    ifelse(x=='1', 0.05,
    ifelse(x=='1L', 0.025,
    ifelse(x=='1U', 0.075,
    ifelse(x=='2', 0.2,
    ifelse(x=='2L', 0.15,
    ifelse(x=='2U', 0.25,
    ifelse(x=='3', 0.4,
    ifelse(x=='3L', 0.35,
    ifelse(x=='3U', 0.45,
    ifelse(x=='4', 0.625,
    ifelse(x=='4L', 0.5625,
    ifelse(x=='4U', 0.6875,
    ifelse(x=='5', 0.875,
    ifelse(x=='5L',0.8125,0.9375)))))))))))))))
}




data =
  data %>%
  dplyr::select(-matches("^AIMS_REEF_NAME$|^REGION$|^NRM_REGION$|^A_SECTOR$|^FRAME$|^CRUISE_CODE$|$MMP_REGION$")) %>%
  ## dplyr::rename(YEAR=one_of("REPORT_YEAR","YEAR"),
  ##               SURVEY_DATE=one_of("SAMPLE_DATE","SURVEY_DATE")) %>%
  filter(REPORT_YEAR>1985) %>%
  mutate(Cover=LTMP_calcPercent(LIVE_CORAL)) %>%
  mutate(
    fYEAR=REPORT_YEAR,
    across(c(TOW_SEQ_NO, fYEAR), function(x) factor(as.character(x)))
    ## DATE=as.Date(SURVEY_DATE, format='%Y-%m-%d')
  ) %>%
  dplyr::select(-matches("^SAMPLE_CLASS$")) %>%
  #group_by(REEF, REPORT_YEAR, fYEAR, DATE) %>%
  #summarise(Cover=mean(Cover, na.rm=TRUE),
  #          Tows=length(unique(TOW_SEQ_NO))) %>%
  mutate(fTOW=factor(paste0(REEF_NAME,TOW_SEQ_NO)),
    nTows=length(unique(TOW_SEQ_NO))) %>%
  ungroup %>%
  ## left_join(data.spatial %>% dplyr::select(-Latitude,-Longitude)) %>%
  mutate(VARIABLE='HC') %>%
  suppressMessages() %>%
  suppressWarnings()

## save(data,  file=paste0(DATA_PATH,'processed/',RDATA_FILE))

m <- vector('list', length(unique(data$REEF_NAME)))
names(m) <- unique(data$REEF_NAME)

for (r in unique(data$REEF_NAME)) {
  dat <- data |> filter(REEF_NAME == r) |> droplevels()
  m[[r]] <- fitModel(dat, REEF_NAME = r)
}
m <- do.call('rbind', m)
save(m, file = "../data/modelled/m.RData")

m1 <- m |> as_tibble() |>
  dplyr::select(posteriors.year) |> 
  unnest(posteriors.year)
## exclude those without fYEAR == "2023"
rfs <- m1 |>
  filter(.draw ==  1) |> 
  mutate(YEAR = as.numeric(as.character(fYEAR))) |>
  group_by(REEF_NAME) |> 
  mutate(COUNT = sum(YEAR > (final_year - 3)),
    maxYear =  max(YEAR)) |> 
  filter(maxYear ==  final_year & COUNT > 1) |>
  droplevels() |>
  pull(REEF_NAME) |>
  unique()
m2 <- m1 |>
  filter(REEF_NAME %in% rfs) |>
  mutate(YEAR = as.numeric(as.character(fYEAR))) |>
  group_by(.draw, REEF_NAME) |>
  arrange(desc(YEAR)) |>
  slice(1:2) |>
  summarise(
    Diff = -diff(value),
    YearDiff = -diff(YEAR),
    AnnualDiff = Diff / YearDiff,
    RefYear = first(YEAR),
    Cover = value[1]
  )

## At some point during 2024, Mike needed the differences between 2024
## and the previous visit for each reef
Diffs_4_Mike <- m2 |>
  ungroup() |> 
  dplyr::select(-YearDiff, -RefYear, -Diff, -Cover) |>
  posterior::as_draws() |>
  group_by(REEF_NAME) |> 
  tidybayes::summarise_draws(Median =  median,
    Mean =  mean,
    HDInterval::hdi,
    `P>0` = ~mean(.x>0),
    `P<0` = ~mean(.x<0))
save(Diffs_4_Mike, file = "../data/modelled/Diffs_4_Mike.RData")

## For the coral change map, we need a data set that has the cover for
## the final_year along with the annual change between this final_year
## and the previous survey (provided it is within the last 3 years).
coral_cover_and_change <- m2 |>
  ungroup() |>
  mutate(REPORT_YEAR = RefYear) |> 
  dplyr::select(-YearDiff, -RefYear, -Diff) |>
  posterior::as_draws() |>
  group_by(REPORT_YEAR, REEF_NAME) |> 
  tidybayes::summarise_draws(Median =  median,
    Mean =  mean,
    HDInterval::hdi,
    `P>0` = ~mean(.x>0),
    `P<0` = ~mean(.x<0))
save(coral_cover_and_change, file = "../data/modelled/coral_cover_and_change.RData")



LTMP_analysis_splits <- function(DOMAIN_CATEGORY, data, VARIABLE){
  data = data %>% filter(VARIABLE==VARIABLE) %>% droplevels
  analysis.splits=NULL
  if (DOMAIN_CATEGORY=='reef') {
    splits <- tribble(
      ## MODEL_GROUPS      - include fYEAR*fGROUP in models
      ## SPLIT_REEF_ZONE   - fit separate models per reef zone
      ## MODEL_REEF_ZONE   - append REEF_NAME with REEF_ZONE (and have this flow through SITE and TRANSECT) so that each ZONE is a separate reef
      ## SPLIT_DEPTH       - fit separate models per depth
      ## MODEL_DEPTH       - include f(fDEPTH, model='iid') in models
      ## SPLIT_SHELF       - fit separate models per shelf
      ~VARIABLE,   ~MODEL_GROUPS, ~SPLIT_REEF_ZONE, ~MODEL_REEF_ZONE, ~SPLIT_DEPTH,  ~REPORT_DEPTH, ~MODEL_DEPTH, ~SPLIT_SHELF,
      VARIABLE,    TRUE,          TRUE,             FALSE,            TRUE,          TRUE,          FALSE,        FALSE,
      )

  } else {
    splits <- tribble(
      ## MODEL_GROUPS      - include fYEAR*fGROUP in models
      ## SPLIT_REEF_ZONE   - fit separate models per reef zone
      ## MODEL_REEF_ZONE   - append REEF_NAME with REEF_ZONE (and have this flow through SITE and TRANSECT) so that each ZONE is a separate reef
      ## SPLIT_DEPTH       - fit separate models per depth
      ## MODEL_DEPTH       - include f(fDEPTH, model='iid') in models
      ## SPLIT_SHELF       - fit separate models per shelf
      ~VARIABLE,   ~MODEL_GROUPS, ~SPLIT_REEF_ZONE, ~MODEL_REEF_ZONE, ~SPLIT_DEPTH,  ~REPORT_DEPTH, ~MODEL_DEPTH, ~SPLIT_SHELF,
      VARIABLE,    FALSE,         FALSE,             TRUE,            FALSE,          TRUE,          TRUE,        TRUE,
      )
  }
  if (data %>% pull(fGROUP) %>% unique %>% length == 1) splits$MODEL_GROUPS=FALSE
  if (data %>% pull(REEF_ZONE) %>% unique %>% length == 1) splits$MODEL_REEF_ZONE=FALSE
  if (data %>% pull(fDEPTH) %>% unique %>% length == 1) splits$MODEL_DEPTH=FALSE
  return(splits)
}



fitModel <- function(data, REEF_NAME) {
  
  model.hierarchy = list(
    GBR=formula(~f(NRM_region, model='iid')+f(REEF, model='iid')+f(REEF_YEAR, model='iid')),
    Zones=formula(~f(REEF, model='iid')+f(REEF_YEAR, model='iid')),
    nrm=formula(~f(REEF, model='iid')+f(REEF_YEAR, model='iid')),
    Bioregions=formula(~f(REEF, model='iid')+f(REEF_YEAR, model='iid')),
    Sectors=formula(~f(REEF, model='iid')+f(REEF_YEAR, model='iid')),
    reef=formula(~1)
  )
  form=model.hierarchy[[DOMAIN_CATEGORY]]
  form=update.formula(form, Cover ~.+fYEAR)


  variables <- c('HC')
  newdata.year=vector('list', length(variables))
  names(newdata.year) <- variables
  newdata.yearcomp <- newdata.year
  posteriors.year <- newdata.year
  if (DOMAIN_CATEGORY=='reef') newdata.yeargroup <- newdata.year
  for (v in variables) {
    cat(paste0('Variable:', v, '\n\n'))
    ## Determine whether the data should be split
    depth_col=c(fGROUP=NA, fDEPTH=NA, REEF_ZONE=NA)  #to enable us to determine whether the fDEPTH column exists and if not what it should include
    data <- data %>%
      add_column(!!!depth_col[!names(depth_col) %in% names(.)])  # add fDEPTH if it does not already exist
    analysis.splits = LTMP_analysis_splits(DOMAIN_CATEGORY=DOMAIN_CATEGORY, data=data, VARIABLE=v)
    if (DOMAIN_CATEGORY == 'Sectors') analysis.splits$SPLIT_SHELF <- FALSE
    data.sub <- data %>%
      filter(VARIABLE==v) %>% droplevels %>%
      mutate(
        INCLUDE_REEF_ZONE=analysis.splits$MODEL_REEF_ZONE,
        INCLUDE_DEPTH=analysis.splits$MODEL_REEF_ZONE,
        SPLIT_REEF_ZONE=analysis.splits$SPLIT_REEF_ZONE,
        SPLIT_DEPTH=analysis.splits$SPLIT_DEPTH,
        SPLIT_SHELF=analysis.splits$SPLIT_SHELF,
        SPLITS = paste0(
               ifelse(SPLIT_REEF_ZONE, as.character(REEF_ZONE),'_'),
               ifelse(SPLIT_DEPTH, as.character(fDEPTH),'_'),
               ifelse(SPLIT_SHELF, as.character(Shelf), '_')),
        ## Prepare the spatial hierarchy
        REEF=factor(ifelse(INCLUDE_REEF_ZONE, paste(REEF, REEF_ZONE),REEF)),
        REEF_YEAR = interaction(REEF, fYEAR)#,
        ## SITE_NO=factor(ifelse(INCLUDE_DEPTH, paste(REEF,SITE_NO,fDEPTH), paste(REEF,SITE_NO))),
        ## TRANSECT_NO=factor(paste(SITE_NO, TRANSECT_NO))
      ) %>%
      dplyr::select(-INCLUDE_REEF_ZONE, -INCLUDE_DEPTH, -SPLIT_REEF_ZONE, -SPLIT_DEPTH, -SPLIT_SHELF)

    ## Determine whether taxonomic groups should be modelled
    if(analysis.splits$MODEL_GROUPS) {
      newdata.yeargroup[[v]] <- newdata.year[[v]]
      form <- update.formula(form, . ~ fYEAR*fGROUP + .)
      fgroups <- TRUE
    } else {fgroups <- FALSE}
    ## Determine whether depth should be modelled as a random effect
    if (analysis.splits$MODEL_DEPTH) form <- update.formula(form, .~. + f(fDEPTH, model='iid'))
    ## Establish lists according to splits
    SPLITS <- unique(data.sub$SPLITS)
    ## In addition to inshore and offshore, there is a desire to have an additional
    ## 'All' SHELF that is both of them together.
    ## That being the case, we will define a new SPLIT (All).
    ## This category does not appear in the data, so we need to be a bit careful when filtering etc
    SHELF_ALL <- 'All'
    if (analysis.splits$SPLIT_SHELF) {
      ss <- unique(gsub('(.*)_(.*)_.*','\\1_\\2_All', SPLITS))
      SPLITS <- c(SPLITS, ss)
    } else ss <- unique(SPLITS)
    newdata.year[[v]] = vector('list', length(SPLITS))
    names(newdata.year[[v]]) <- SPLITS
    newdata.yearcomp[[v]] <- newdata.year[[v]]

    ## if (length(unique(data$fYEAR))<2) form=update.formula(form, .~.-fYEAR - fYEAR:fGROUP)
    cat(paste0('Fitting with the formula: ', Reduce(paste, deparse(form)), '\n'))
    for (d in SPLITS) {
    ##   if (d %in% ss) {  #provide a mechanism for including all shelf positions
    ##     data.group <- data.sub %>%
    ##       mutate(SPLITS = gsub(paste0('(',unique(data.sub$Shelf),')', collapse='|'),SHELF_ALL, SPLITS, perl=TRUE),
    ##         Shelf = 'All') %>%
    ##       filter(SPLITS == d) %>%
    ##       droplevels
    ##   } else {
        data.group <- data.sub %>% filter(SPLITS==d) %>% droplevels()
      ## }
      ## if (length(unique(data.group$fYEAR))<2) form=remove_terms(form, 'fYEAR') #form=update.formula(form, .~.-fYEAR - fYEAR:fGROUP)
      if (length(unique(data.group$fYEAR))<2) {
        return(list(newdata.year =  NULL, yeardata.yearcomp =  NULL, posteriors =  NULL))
      }
      mod_plus = fitMantaModel(data.group, form)
      ## zone=data.group %>% pull(REEF_ZONE) %>% unique  #gsub('(.*)\\.[0-9]','\\1',d)
      ## zone=ifelse(length(zone)==1, zone, 'Various')
      ## zone=NA
      ## depth=data.group %>% pull(fDEPTH) %>% unique() %>% as.character() %>% as.numeric() %>% mean(na.rm=TRUE) #mean(as.numeric(gsub('.*\\.([0-9])','\\1',d)), na.rm=TRUE)
      ## depth=NA
      zone=data.group %>% pull(REEF_ZONE) %>% unique  #gsub('(.*)\\.[0-9]','\\1',d)
      zone=ifelse(length(zone)==1, zone, 'Various')
      depth=data.group %>% pull(fDEPTH) %>% unique() %>% as.character() %>% as.numeric() %>% mean(na.rm=TRUE) #mean(as.numeric(gsub('.*\\.([0-9])','\\1',d)), na.rm=TRUE)
      depth=ifelse(is.nan(depth),NA, depth)
      ## shelf=data.group %>% pull(Shelf) %>% unique
    shelf <- NA
      shelf=ifelse(length(shelf)==1 & !grepl('.*All', shelf), as.character(shelf), 'Various')
      newdata.year[[v]][[d]] = data.frame(VARIABLE=v, DEPTH=depth, REEF_ZONE=zone, SHELF=shelf, mod_plus[['newdata.year']]) %>%
        as_tibble()
      newdata.yearcomp[[v]][[d]] = data.frame(VARIABLE=v, DEPTH=depth, REEF_ZONE=zone, SHELF=shelf, mod_plus[['newdata.yearcomp']]) %>%
        as_tibble()
      posteriors.year[[v]][[d]] <- data.frame(VARIABLE = v, DEPTH = depth, REEF_ZONE = zone, SHELF = shelf, REEF_NAME =  REEF_NAME,
        mod_plus[['posteriors']]) |> as_tibble()
    }
    newdata.year[[v]] = do.call('rbind', newdata.year[[v]])
    newdata.yearcomp[[v]] = do.call('rbind', newdata.yearcomp[[v]])
    posteriors.year[[v]] = do.call('rbind', posteriors.year[[v]])
  }
  newdata.year = do.call('rbind', newdata.year)
  newdata.yearcomp = do.call('rbind', newdata.yearcomp)
  posteriors.year = do.call('rbind', posteriors.year)

  list(newdata.year =  newdata.year, newdata.yearcomp =  newdata.yearcomp, posteriors.year =  posteriors.year)
}
  

fitMantaModel <- function(data.group, form){
## print('Before')
  ## print(levels(data.group$fYEAR))
  data.group = data.group %>% arrange(fYEAR) %>% mutate(fYEAR=fct_relevel(fYEAR, rev))
## print('After')
## print(levels(data.group$fYEAR))
  ## Start by appending prediction data with the observed data
  newdata = with(data.group, data.frame(fYEAR=levels(fYEAR),
                                        REEF=NA, Cover=NA, Tows=NA, fTOW=NA))
  ## print(newdata)
  ## Add back on the mean survey date so that the modelled values can be plotted on an x-axis reflecting sampling time
  ## Will need to think about what to do about higher spatial aggregates (domains)
  newdata = newdata %>%
    mutate(REPORT_YEAR=as.numeric(as.character(fYEAR)))##  %>%
    ## left_join(data.group %>%
    ##           group_by(fYEAR) %>%
    ##           summarise(DATE=if_else(length(unique(REEF))==1, mean(DATE, na.rm=TRUE),  as.Date(paste0(unique(fYEAR), '-01-01')))) %>%
    ##           dplyr::select(fYEAR, DATE) %>%
    ##           distinct)
  newdata1 <- newdata
  ## print(newdata)
  ## print(levels(newdata$fYEAR))
  data.pred = data.group %>% bind_rows(newdata) %>%
    mutate(Cover=ifelse(Cover==0,0.001,ifelse(Cover==1, 0.999,Cover)),
           fYEAR=fct_relevel(fYEAR,levels(data.group$fYEAR)))
  ## print(levels(data.pred$fYEAR))
  ## Fit the model
  ##inla.setOption("enable.inla.argument.weights", TRUE)
  mod.inla = inla(form,
                  data=data.pred, family='beta',
                  #weights=data.pred$Tows,
                  control.family=list(link='logit'),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list(config = TRUE, dic=TRUE, cpo=TRUE, waic=TRUE)
                  )
  newdata <- cbind(newdata,
                   mod.inla$summary.fitted.values[(nrow(data.group)+1):nrow(data.pred),]) %>%
    rename(lower=`0.025quant`, upper=`0.975quant`, median=`0.5quant`)
  ##if you want to express year diffs as Risk ratio (instead of odds ratio):
  ##RR = OR/((1-Pr)+(Pr*OR)), where Pr is prob of second year
  RR <- function(OR, Pr) {
    OR/((1-Pr)+(Pr*OR))
  }

  newdata.yearcomp <- data.pred %>% tidyr::expand(fYEAR) %>%
    mutate(fYEAR=paste0(first(fYEAR),'-',fYEAR)) %>%
    slice(-1) %>%
    bind_cols(
      mod.inla$summary.fixed[-1,] %>% mutate(across(where(is.numeric), function(x) 1/exp(x)))
    ) %>%
    dplyr::rename(YearComp=fYEAR, lower=`0.975quant`, upper=`0.025quant`, median=`0.5quant`) %>%
    dplyr::select(YearComp, mean, sd, lower, median, upper, mode, kld)

  ## full posteriors per year
  draws <- inla.posterior.sample(n=1000, mod.inla, seed=123)  
  cellmeans <- sapply(draws, function(x)
    x[[2]][(nrow(data.pred)-nrow(newdata1)+1):nrow(data.pred)]) 
  posteriors <- newdata1 %>%
    dplyr::select(fYEAR) %>%
    cbind(plogis(cellmeans)) %>%
    pivot_longer(cols = matches('[0-9]'), names_to = 'Rep') %>%
    mutate(.draw = as.integer(Rep)) %>%
    dplyr::select(-Rep) %>%
    ## left_join(obs_data %>%
    ##           dplyr::select(fYEAR) %>%
    ##           distinct()) %>% 
    suppressWarnings() %>%
    suppressMessages()
  list(mod=mod.inla, newdata.year=newdata, newdata.yearcomp=newdata.yearcomp, posteriors =  posteriors)
}
