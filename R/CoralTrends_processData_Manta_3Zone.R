cat(paste0("\nProcessing annual report ",
  domain, " data\n====================================\n"),
  append = TRUE)

library(sf)
source('CoralTrends_functions.R')
## CoralTrends_checkPackages()

cat(paste0("\t- adding pcode modifications\n"),
  append = TRUE)
## ---- processManta_getData
pcode_mod <- read.csv('../parameters/P_CODE_MOD.csv', strip.white=TRUE)
manta <- readRDS(file = paste(data_path, "data_", domain, ".rds"))
manta <- manta |>
  left_join(pcode_mod, by = "REEF_NAME")
## ----end

#####################################################################
## Data are collected per Manta Tow.                               ##
## The most appropriate unit for these analyses is the reef level. ##
## - Only include data collected after 1985                        ##
## - Convert data into percent cover                               ##
## - Summarize the data to the reef/year level                     ##
#####################################################################
cat(paste0("\t- prepare spatial components of data\n"),
  append = TRUE)
## ---- processManta_towlevel
manta_tow <- manta |>
  filter(REPORT_YEAR > 1985) |>
  mutate(Cover = CoralTrends_calcPercent(LIVE_CORAL)) |>
  mutate(Latitude = LATITUDE, Longitude = LONGITUDE) |>
  CoralTrends_calc3ZoneLocation() |>
  filter(!is.na(Region)) |>
  droplevels() |>
  mutate(Zone = Region) |>
  ## mutate(Region = factor(Region, levels = c('Northern GBR', 'Central GBR', 'Southern GBR')),
  ##        Zone=Region) |>
  mutate(Year = factor(REPORT_YEAR)) |>
  ## mutate(Cover=ifelse(Cover==0,0.0001,Cover)) |>
  mutate(Cover = ifelse(Cover == 0,0.0001, Cover)) |>
  mutate(P_CODE.mod = factor(ifelse(is.na(P_CODE.mod), 'Other', P_CODE.mod))) |>
  as.data.frame()
## ----end

## ---- processManta_aggregateToReef
manta_sum <- manta |>
  filter(REPORT_YEAR > 1985) |>
  mutate(Cover = CoralTrends_calcPercent(LIVE_CORAL)) |>
  group_by(P_CODE.mod,
    ## A_SECTOR,
    ## SHELF,
    REEF_NAME,
    AIMS_REEF_NAME,
    ## REEF_ID,
    REPORT_YEAR) |>
  summarise(Cover = mean(Cover, na.rm = TRUE),
    CoverCat.median = median_cover_cat(LIVE_CORAL),
    Tows = length(unique(TOW_SEQ_NO)),
    Latitude = mean(LATITUDE, na.rm = TRUE),
    Longitude=mean(LONGITUDE, na.rm = TRUE)) |>
  mutate(Cover_from_cat = CoralTrends_calcPercent(CoverCat.median)) |>
  ungroup() |>
  group_by(REEF_NAME) |>
  ungroup() |>
  CoralTrends_calc3ZoneLocation() |>
  filter(!is.na(Region)) |>
  droplevels() |>
  mutate(Zone = Region) |>
  ## mutate(Region = factor(Region, levels = c('Northern GBR', 'Central GBR', 'Southern GBR')),
  ##   Zone = Region) |>
  as.data.frame()
## ----end

## ---- Compare tow and reef level
## bfun <- function(x) {
##     if (length(x$Cover)>3 & !all(x$Cover==x$Cover[1])) {
##         betareg::betareg(Cover ~ 1, data=x) %>% coef() %>% `[[`(1) %>% plogis()
##     } else {
##         x$Cover[1]
##     }
## }
## manta_stats <- manta_tow |>
##     mutate(RN = REEF_NAME, YR = Year) |>
##     group_by(REEF_NAME, Year) |>
##     nest() |>
##     mutate(Mean = map_dbl(data, ~mean(.x$Cover)),
##            Median = map_dbl(data, ~median(.x$Cover)),
##            Logis = map_dbl(data, ~.x$Cover |> gtools::logit() |> mean() |> plogis()),
##            Beta = map_dbl(data, ~bfun(.x))
##            ) |>
##     dplyr::select(-data) |>
##     unnest(cols = c()) |>
##     relocate(Mean, Beta, .after = last_col()) |>
##     arrange(Year) |>
##     ungroup() |>
##     as.data.frame()
## ----end

cat(paste0("\t- save the processed data\n"),
  append = TRUE)
## ---- Save processed data
cat(paste0("\t\t- NumberOfReefs=",length(unique(manta_sum$REEF_NAME)),
                  "\n\t\t- NumberOfSurveys=",nrow(manta_sum)), append = TRUE)
saveRDS(manta_sum, file = paste0(data_path, "manta_sum_", domain, ".rds"))
saveRDS(manta_tow, file = paste0(data_path, "manta_tow_", domain, ".rds"))
## ---- end

cat(paste0("\nProcessing of ",
  domain, " data complete!\n"),
  append = TRUE)
