source('CoralTrends_functions.R')
CoralTrends_checkPackages()
source('CoralTrends_config.R')

## ---- load Data
load('../data/processed/manta.sum.RData')
## load(file='../data/modelled/dat.gbr.RData')
load(file='../data/modelled/mod.northern_brms.beta.ry.disp.RData')
load(file='../data/modelled/mod.central_brms.beta.ry.disp.RData')
load(file='../data/modelled/mod.southern_brms.beta.ry.disp.RData')
dat.gbr <- rbind(dat.northern_brms.beta.ry.disp,
                 dat.central_brms.beta.ry.disp,
                 dat.southern_brms.beta.ry.disp
                 )
#load(file='../data/modelled/last_year.RData')
## load(file='../data/modelled/mod.gbr.RData')
mod.gbr <- NULL

load(file='../data/spatial/whagbr.RData')
load(file='../data/spatial/whagbr.n.RData')
load(file='../data/spatial/whagbr.c.RData')
load(file='../data/spatial/whagbr.s.RData')
load(file='../data/spatial/qld.RData')

## ----end

## In 2023, Mike needs specific coral recovery info:

## - Northern 2017 - 2022 
## - Central 2012 - 2016, 2019 - 2022
## - Southern 2011 - 2017, 2011 - 2021

library(tidybayes)
## Northern
emmeans.northern_brms.beta.ry.disp = emmeans(mod.northern_brms.beta.ry.disp,
                                             ~Year,
                                             at = list(Year = c(2017, 2022)),
                                             type='link') %>%
    gather_emmeans_draws()
emmeans.northern_brms.beta.ry.disp %>% head() 
northern.change <- emmeans.northern_brms.beta.ry.disp %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2022` - `2017`,
           Changep = Change/`2017`,
           ChangeP = exp((log(`2022`) - log(`2017`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2022 - 2017),
           AnnualChangePP = exp((log(`2022`) - log(`2017`))/(2022 - 2017)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2017`, -`2022`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    ungroup() %>%
    group_by(name) %>%
    median_hdci(value) %>%
    as.data.frame()
northern.change 
save(northern.change, file = paste0("../data/modelled/northern.change.RData"))


## Central
emmeans.central_brms.beta.ry.disp = emmeans(mod.central_brms.beta.ry.disp,
                                             ~Year,
                                             at = list(Year = c(2012, 2016)),
                                             type='link') %>%
    gather_emmeans_draws()
emmeans.central_brms.beta.ry.disp %>% head() 
central.change <- emmeans.central_brms.beta.ry.disp %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2016` - `2012`,
           Changep = Change/`2012`,
           ChangeP = exp((log(`2016`) - log(`2012`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2016 - 2012),
           AnnualChangePP = exp((log(`2016`) - log(`2012`))/(2016 - 2012)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2012`, -`2016`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    ungroup() %>%
    group_by(name) %>%
    median_hdci(value) %>%
    as.data.frame()
central.change 

emmeans.central_brms.beta.ry.disp1 = emmeans(mod.central_brms.beta.ry.disp,
                                             ~Year,
                                             at = list(Year = c(2019, 2022)),
                                             type='link') %>%
    gather_emmeans_draws()
emmeans.central_brms.beta.ry.disp1 %>% head() 
central.change1 <- emmeans.central_brms.beta.ry.disp1 %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2022` - `2019`,
           Changep = Change/`2019`,
           ChangeP = exp((log(`2022`) - log(`2019`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2022 - 2019),
           AnnualChangePP = exp((log(`2022`) - log(`2019`))/(2022 - 2019)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2019`, -`2022`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    ungroup() %>%
    group_by(name) %>%
    median_hdci(value) %>%
    as.data.frame()
central.change1 

save(central.change, central.change1,
     file = paste0("../data/modelled/northern.change.RData"))

## Southern
emmeans.southern_brms.beta.ry.disp = emmeans(mod.southern_brms.beta.ry.disp,
                                             ~Year,
                                             at = list(Year = c(2011, 2017)),
                                             type='link') %>%
    gather_emmeans_draws()
emmeans.southern_brms.beta.ry.disp %>% head() 
southern.change <- emmeans.southern_brms.beta.ry.disp %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2017` - `2011`,
           Changep = Change/`2011`,
           ChangeP = exp((log(`2017`) - log(`2011`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2017 - 2011),
           AnnualChangePP = exp((log(`2017`) - log(`2011`))/(2017 - 2011)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2011`, -`2017`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    ungroup() %>%
    group_by(name) %>%
    median_hdci(value) %>%
    as.data.frame()
southern.change 

emmeans.southern_brms.beta.ry.disp1 = emmeans(mod.southern_brms.beta.ry.disp,
                                             ~Year,
                                             at = list(Year = c(2011, 2021)),
                                             type='link') %>%
    gather_emmeans_draws()
emmeans.southern_brms.beta.ry.disp1 %>% head() 
southern.change1 <- emmeans.southern_brms.beta.ry.disp1 %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2021` - `2011`,
           Changep = Change/`2011`,
           ChangeP = exp((log(`2021`) - log(`2011`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2021 - 2011),
           AnnualChangePP = exp((log(`2021`) - log(`2011`))/(2021 - 2011)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2011`, -`2021`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    ungroup() %>%
    group_by(name) %>%
    median_hdci(value) %>%
    as.data.frame()
southern.change1 

save(southern.change, southern.change1,
     file = paste0("../data/modelled/northern.change.RData"))


## Comparisons

change.northern <- emmeans.northern_brms.beta.ry.disp %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2022` - `2017`,
           Changep = Change/`2017`,
           ChangeP = exp((log(`2022`) - log(`2017`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2022 - 2017),
           AnnualChangePP = exp((log(`2022`) - log(`2017`))/(2022 - 2017)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2017`, -`2022`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    mutate(Zone = 'Northern')


change.central <- emmeans.central_brms.beta.ry.disp %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2016` - `2012`,
           Changep = Change/`2012`,
           ChangeP = exp((log(`2016`) - log(`2012`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2016 - 2012),
           AnnualChangePP = exp((log(`2016`) - log(`2012`))/(2016 - 2012)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2012`, -`2016`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    mutate(Zone = 'Central 1')


change.central1 <- emmeans.central_brms.beta.ry.disp1 %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2022` - `2019`,
           Changep = Change/`2019`,
           ChangeP = exp((log(`2022`) - log(`2019`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2022 - 2019),
           AnnualChangePP = exp((log(`2022`) - log(`2019`))/(2022 - 2019)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2019`, -`2022`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    mutate(Zone = 'Central 2')

change.southern <- emmeans.southern_brms.beta.ry.disp %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2017` - `2011`,
           Changep = Change/`2011`,
           ChangeP = exp((log(`2017`) - log(`2011`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2017 - 2011),
           AnnualChangePP = exp((log(`2017`) - log(`2011`))/(2017 - 2011)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2011`, -`2017`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    mutate(Zone = 'Southern 1')

change.southern1 <- emmeans.southern_brms.beta.ry.disp1 %>%
    mutate(.value = plogis(.value)) %>%
    pivot_wider(id_cols = everything(),
                names_from = "Year",
                values_from = .value) %>%
    mutate(Change = `2021` - `2011`,
           Changep = Change/`2011`,
           ChangeP = exp((log(`2021`) - log(`2011`))),
           ChangePP = 100*(ChangeP - 1),
           AnnualChange = 100 * Change/(2021 - 2011),
           AnnualChangePP = exp((log(`2021`) - log(`2011`))/(2021 - 2011)),
           AnnualChangePP = 100*(AnnualChangePP),
           Change = 100 * Change) %>%
    dplyr::select(-`2011`, -`2021`, -ChangeP, -Changep) %>%
    pivot_longer(cols = c(Change, ChangePP, AnnualChange, AnnualChangePP)) %>%
    mutate(Zone = 'Southern 2')

Xmat <- cbind("N vs C1" = c(1,-1,0,0,0),
              "N vs C2" = c(1,0,-1,0,0),
              "N vs S1" = c(1,0,0,-1,0),
              "N vs S2" = c(1,0,0,0,-1),
              "C1 vs C2" = c(0,1,-1,0,0),
              "C1 vs S1" = c(0,1,0,-1,0),
              "C1 vs S2" = c(0,1,0,0,-1),
              "C2 vs S1" = c(0,0,1,-1,0),
              "C2 vs C2" = c(0,0,1,0,-1),
              "S1 vs S2" = c(0,0,0,1,-1))

change.all <- bind_rows(change.northern,
                        change.central,
                        change.central1,
                        change.southern,
                        change.southern1) %>%
    filter(name %in% c("AnnualChange", "AnnualChangePP")) %>%
    group_by(.iteration, .draw, name) %>%
    nest() %>%
    mutate(A = map(.x = data,
                   .f = ~ as.data.frame(.x %>% pull(value) %*% Xmat) %>%
                       pivot_longer(cols = everything(), names_to = "Comparison")
                   )
           ) %>%
    unnest(A) %>%
    ungroup() %>%
    group_by(name, Comparison) %>%
    median_hdci(value)

    summarise(A = value %*% Xmat)
change.all %>% dplyr::select(-A)
