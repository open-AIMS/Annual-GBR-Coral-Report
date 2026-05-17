cat(paste0("\nFit the model for the ",
  domain, " data\n====================================\n"),
  append = TRUE)

source('CoralTrends_functions.R')
source('CoralTrends_model_functions.R')
## CoralTrends_checkPackages()
## source('CoralTrends_config.R')
library(INLA)
library(emmeans)
## library(DHARMa)
library(brms)


## Need to convert this into a nested data frame to be consistent
## with the other LTMP dashboard model outputs
data <- tibble(
  VARIABLE = "HC",
  model_type = "Cover",
  family_type = "beta",
  model_response = "Cover",
  sub_model = NA,
  lab = paste0("annual_report_region_", domain),
  analysis_splits = NA,
  splits = NA,
  newdata =  NA,
  raw_sum = NA,
  label = paste0("annual_report_region_", domain,
    "_HC_beta_ _ _ _Cover_ "),
  fig_label = paste0("annual_report_region_", domain,
    "_HC_ _ _ _Cover_ ")
  )


cat(paste0("1. get the data\n"), #-------------------------------------
  append = TRUE)
manta_sum <- readRDS(file = paste0(data_path, "manta_sum_", domain, ".rds"))
manta_tow <- readRDS(file = paste0(data_path, "manta_tow_", domain, ".rds"))


cat(paste0("2. select the model type\n"), #-------------------------------------
  append = TRUE)
## Define the specific model type.
## Note.  In the code that follows, only one model type "BRMS beta ry disp"
## will be defined.  For the others, refer to ../scripts/CoralTrends_stan_Manta_3Zone.R
## This is the old codebase for this project
models <- c(
    ##'BRMS beta vanilla',
    #'BRMS beta disp'#,
    'BRMS beta ry disp'#,
    ##'glmmTMB beta vanilla',
    #'glmmTMB_tow beta disp',
    #'glmmTMB_tow beta ry disp'
    ##'BRMS ordinal'
)
zone <- domain  ## note, this use to just house 'northern, central or southern'

#######################################################################
## Tow level data                                                    ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _raw_data.rds")           ##
#######################################################################
cat(paste0("3. prepare the tow level data\n"),
  append = TRUE)
data <- data |> CoralTrends_prepare_tow_level_data(manta_tow)

#######################################################################
## Define model formula                                              ##
#######################################################################
cat(paste0("4. define the model formula\n"), # -----------------------
  append = TRUE)
data <- data |> CoralTrends_model_formula()

#######################################################################
## Define model priors                                               ##
#######################################################################
cat(paste0("5. define the model priors\n"), #-------------------------
  append = TRUE)
data <- data |> CoralTrends_model_priors()

#######################################################################
## Fit the model                                                     ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _.rds")                   ##
#######################################################################
cat(paste0("6. fit the model\n"), #-----------------------------------
  append = TRUE)
data <- data |> CoralTrends_fit_model()

## mod <- readRDS(file = paste0(data_path, "modelled/annual_report_region_", domain,
##   "_HC_beta_ _ _ _Cover_ .rds"))



data <- data |> CoralTrends_prepare_for_contrasts()


data <- data |> mutate(selected = TRUE)

#######################################################################
## Get cell means posteriors                                         ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _year_posteriors.rds")    ##
#######################################################################
cat(paste0("7. get cell means posteriors\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_cellmeans_posteriors()

#######################################################################
## Get cell means summaries                                          ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _year_sum.rds")           ##
#######################################################################
cat(paste0("8. get cell means summaries\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_cellmeans_sum()

#######################################################################
## Get contrast posteriors                                           ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _yearcomp_posteriors.rds")##
#######################################################################
cat(paste0("7. get contrast posteriors\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_contrast_posteriors()

#######################################################################
## Get contrast summaries                                            ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _yearcomp_sum.rds")       ##
#######################################################################
cat(paste0("8. get contrast summaries\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_contrast_sum()

#######################################################################
## Get all contrast posteriors                                       ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _all_yearcomp_posteriors.rds")##
#######################################################################
cat(paste0("9. get all contrast posteriors\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_all_contrast_posteriors()

#######################################################################
## Get contrast summaries                                            ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _all_yearcomp_sum.rds")       ##
#######################################################################
cat(paste0("10. get all contrast summaries\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_all_contrast_sum()


#######################################################################
## Generate raw summaries                                            ##
## - data are stored in data_path, "modelled/annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ _raw_sums.rds")           ##
#######################################################################
cat(paste0("12. generate raw summaries\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_raw_summary()


#######################################################################
## Generate raw summary plots                                        ##
## - data are stored in data_path, "figures/gg_raw_sum_annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ .rds")                    ##
#######################################################################
cat(paste0("13. generate raw summaries plot\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_raw_summary_plots()

#######################################################################
## Generate raw summary plots                                        ##
## - data are stored in data_path, "figures/gg_annual_report_region_", ##
##       domain, "_COVER_beta_ _ _ _Cover_ .rds")                    ##
#######################################################################
cat(paste0("14. generate raw summaries plot\n"), #------------------------
  append = TRUE)
data <- data |> CoralTrends_partial_plots()

saveRDS(data, file = paste0(data_path, "modelled/annual_report_region_", domain,
  ".rds"))

cat(paste0("Complete!\n"), #------------------------
  append = TRUE)

## data <- readRDS(file = paste0(data_path, "modelled/annual_report_region_", domain,
##   ".rds"))

## raw_sum <- readRDS(file = paste0(data_path, "modelled/annual_report_region_", domain,
##   "_HC_beta_ _ _ _Cover_ _raw_sum.rds"))

year_sum <- readRDS(file = paste0(data_path, "modelled/annual_report_region_", domain,
   "_HC_beta_ _ _ _Cover_ _year_sum.rds"))
