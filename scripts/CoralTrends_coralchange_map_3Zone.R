library(ggsn) #map features
source('CoralTrends_functions.R')
CoralTrends_checkPackages()
source('CoralTrends_config.R')

## Remove this
## final_year <- 2024


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
qld_sf <- st_as_sf(qld)
gbr_sf <- st_read("../data/spatial/Features/Great_Barrier_Reef_Features.shp")
Zones <- c('Northern','Central','Southern')

trafficLightColors = colorRampPalette(c('#ED1C24','#F47721','#F0C918','#B0D235','#00734D'))

## Supplied by Mike
manta.mod <- read_csv("../data/modelled/data-manta-2024.csv") 
## The following could be used as a lookup of lat/long incase the conversion I am using below is no good!
load("../data/primary/reef.lookup.RData")
## manta.mod <- manta.mod |>
##   left_join(
##     reef.lookup |>
##       ungroup() |> 
##       filter(P_CODE %in% c("RM", "RMRAP", "RAP")) |> 
##       dplyr::select(AIMS_REEF_NAME, REEF_LAT, REEF_LONG) |>
##       distinct(),
##     by = c("aims_reef_name" = "AIMS_REEF_NAME")
##   )

if (final_year == 2024) {
    cotsAndBleaching <- read_csv("../data/primary/202324_LTMP_Annual_report_bleaching_COTS 1.csv") %>%
        mutate(Ave_perc_bleach_cat = as.character(ltmp_bleach_cat),
               ## Ave_perc_bleach_cat = ifelse(Ave_perc_bleach_cat == "0", "0+", Ave_perc_bleach_cat),
               ## `in-water_ sample_date` = in.water_.sample_date)
               Aerial_bleaching = as.character(aerial_bleach),
               Aerial_date = as.Date(aerial_date, format = "%d/%m/%Y"))
}
library(magick)
library(png)
water_mark = magick::image_read(path='../parameters/ECM_1280725_v1_AIMS Logo - stacked.jpg') 
## ----end

ggplot() +
  geom_sf(data = gbr_sf) +
  geom_point(data = manta.mod, aes(y = REEF_LAT, x = REEF_LONG))




## Before looking at modelled data, lets make a plot of raw data for the finalyear

## ---- Prepare axes labels and color scales
ewbrks <- seq(144,152,by=2)
nsbrks <- seq(-24,-10,by=2)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste0(x, "°W"), ifelse(x > 0, paste0(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), "°S"), ifelse(x > 0, paste0(x, "°N"),x))))
trafficLightColors = colorRampPalette(c('#ED1C24','#F47721','#F0C918','#B0D235','#00734D'))
## ----end


## Generate a map with multiple overlays
make_plot_2024 <- function() {
  ## Get the data provided by Mike
  ## manta.mod <- get(load("../data/modelled/ltm-modelled-data-from_dashboard_2021_2023.RData")) |>
  manta.mod <- manta.mod |> process_manta_mod_data()
  management_sf <- get_gbr_sf()
  management_df <- get_gbr_df(management_sf)
  tops <- get_gbr_tops(management_df)
  coefs <- get_management_coefs(management_df)
  lgnd.dat <- get_lgnd_pos(coefs)
  ## 1. Base plot
  survey_dates <- cotsAndBleaching |>
    mutate(Date = as.Date(`in-water_sample_date`, format = "%d/%m/%Y")) |> 
    mutate(Region = factor(Region, levels = c("Northern", "Central", "Southern"))) |>
    arrange(Region) |> 
    group_by(Region) |>
    summarise(min = min(Date, na.rm = TRUE), max = max(Date, na.rm = TRUE))

  survey_dates <-
    survey_dates |>
    mutate(survey_date = ifelse(month(min) == month(max),
      paste0(month(min, label = TRUE), " ", year(min)),
      paste0(
        month(min, label = TRUE), " ", year(min), " - ",
        month(max, label = TRUE), " ", year(max)
      )
    )) |>
    dplyr::select(Region, survey_date)
  SurveyDates <- paste0("Suvery dates\n(", survey_dates |> pull(survey_date), ")")
  ## SurveyDates <-  c('Survey dates\n(Oct-Dec 2022)',
  ##                   'Survey dates\n(Oct 2022 - May 2023)',
  ##                   'Survey dates\n(Aug 2022 - May 2023)')
  g.base <- make_base_plot(
    manta.mod, gbr_sf, qld, management_sf,
    tops, lgnd.dat, towns12, Zones, SurveyDates,
    extra_width = 23
  )
  ## 2. Coral change
  g.change <- make_change_plot(manta.mod, management_sf, tops)
  g <- add_change_to_plot(g = g.base, g.change, management_sf, shiftacross = 5.6)
  ## 3. COTS
  cots <- process_cots(cotsAndBleaching)
  g.cots <- make_cots_plot(cots, coefs, tops)
  g <- add_change_to_plot(g = g, g.cots, management_sf, shiftacross = 11.2)
  ## 4a. Bleaching
  bleaching <- process_bleaching(cotsAndBleaching, cots)
  g.bleaching <- make_bleaching_plot(bleaching, coefs, tops, management_sf)
  g <- add_change_to_plot(g = g, g.bleaching, management_sf, shiftacross = 16.8)
  ## 4b. Aerial Bleaching
  bleaching.aerial <- process_aerial_bleaching(cotsAndBleaching)
  g.bleaching.aerial <- make_aerial_bleaching_plot(bleaching.aerial, coefs, tops, management_sf) 
  g <- add_change_to_plot(g = g, g.bleaching.aerial, management_sf, shiftacross = 22.4)
  ## 5. Add the Australia inset
  g.oz <- make_oz_inset(qld_sf)
  g <- add_oz_to_plot(g, g.oz, management_sf, extra_width = 23, inset_width = 5, xmin = 170, ymin = -15)
  ## 6. Add the scalebar and north arrow
  g <- add_scale_and_arrow_to_plot(g)
  ## 7. AIMS watermark
  g <- add_watermark_to_plot(g, water_mark, management_sf, extra_width = 23, oz_width = 5, inset_width = 3)
  g
}
g <- make_plot_2024()

## bb.ratio <- (bb[3]-bb[1])/(bb[4]-bb[2])
ggsave(filename=paste0('../output/figures/FourPlots_2024.pdf'), g, width=15,height=15/1.7)
ggsave(filename=paste0('../output/figures/FourPlots_2024.png'), g, width=15,height=15/1.7, dpi = 600)
## ggsave(filename=paste0('../output/figures/FourPlots2.png'), g1, width=10,height=6.8, dpi=300)


process_manta_mod_data <- function(manta.mod) {
  manta.mod |> 
    mutate(Latitude = -lat_deg - lat_min/60,
      Longitude = long_deg + long_min/60,
      REPORT_YEAR = report_year,
      REEF_NAME = aims_reef_name,
      Cover = median)
}

get_gbr_sf <- function() {
  rbind(
    st_as_sf(whagbr.n),
    st_as_sf(whagbr.c),
    st_as_sf(whagbr.s)
  )
}

get_gbr_df <- function(management_sf) {
  management_sf |>
    st_coordinates() |>
    as_tibble() |> 
    rename(long = X, lat = Y) 
}

get_gbr_tops <- function(management_df) {
  ## Postion of labels above GBR polygons
  management_df |>
    filter(lat > max(lat) - 0.1) |>
    summarize(min = min(long), max = max(long), lat = max(lat))
}

ff <- function(strt,coefs,length) {
  x=strt[1]
  y=strt[2]
  dat <- data.frame(x = x, y = y)
  for (i in 2:length) { 
    y <- y - 0.7
    x <- coefs[1] + y*coefs[2]
    dat <- rbind(dat, data.frame(x = x[[1]], y = y[[1]]))
  }
  dat 
}

get_management_coefs <- function(management_df) {
  dat <- management_df |>
    filter(lat < -15 & lat > -16, long > 145.5)
  coefs <- coef(lm(long ~ lat, data = dat))
  coefs
}
get_lgnd_pos <- function(coefs) {
  strt = c(145.9,-14.5)

  lgnd.dat <- ff(c(coefs[1] + (strt[2]*coefs[2]), strt[2]),coefs,length=5)
  lgnd.dat$lab <- rev(c('>0% - 10%','>10% - 30%','>30% - 50%','>50% - 75%','>75% - 100%'))
  lgnd.dat$Cat <- factor(5:1) 
  lgnd.dat
}


make_base_plot <- function(manta.mod, gbr_sf, qld, management_sf, tops,
                           lgnd.dat, towns12, Zones, SurveyDates,
                           extra_width = 15.5) {
  ## define a quasi grid
  grd.y = data.frame(x=c(140,140,140,140,140,140,140),
    xend=c(146,146,148,150,152,154,154),
    y=c(-11,-13,-15,-17,-19,-21,-23),
    yend=c(-11,-13,-15,-17,-19,-21,-23))
  grd.x = data.frame(x=c(142,144,146,148,150,152,154),
    xend=c(142,144,146,148,150,152,154),
    y=c(-9,-9,-11,-15,-17,-19,-21),
    yend=c(-30,-30,-30,-30,-30,-30,-30))

  ## ---- base plot
  g.base <-
    manta.mod |>
    filter(REPORT_YEAR == final_year) |>
    mutate(cCover = cut(Cover, breaks = c(0, 0.10, 0.30, 0.50, 0.75, 1), labels = 1:5)) |>
    ggplot(, aes(y = Latitude, x = Longitude)) +
    geom_sf(data = gbr_sf |> filter(FEAT_NAME == "Reef"), fill = "grey", colour = "grey70") +
    geom_sf(data = qld_sf, fill = "grey", colour = "grey70") +
    geom_segment(data = grd.y, aes(x = x, y = y, xend = xend, yend = yend),
      size = 0.1, color = "grey80") +
    geom_segment(data = grd.x, aes(x = x, y = y, xend = xend, yend = yend),
      size = 0.1, color = "grey80") +
    annotate(geom = "text", x = -Inf, y = seq(-23, -11, by = 2),
      label = degreeLabelsNS(seq(-23, -11, by = 2)), hjust = -0.1, parse = TRUE, size = 3) +
    annotate(geom = "text", y = -Inf, x = seq(142, 154, by = 2),
      label = degreeLabelsEW(seq(142, 154, by = 2)), vjust = -0.3, parse = TRUE, size = 3) +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.5) +
    geom_point(aes(y=Latitude,x=Longitude,fill=cCover),
      size=3, alpha=1,shape=21,color='black', show.legend = TRUE)  +
    scale_fill_manual('Hard coral cover', breaks=1:5,
      values=trafficLightColors(5),
      labels=rev(c('>0% - 10%','>10% - 30%','>30% - 50%','>50% - 75%','>75% - 100%')),
      limits=factor(1:5)) +
    annotate(geom='rect', xmin=tops$min, xmax=tops$max,
      ymin=tops$lat+0.01, ymax=tops$lat+0.7, fill='white', alpha=0.8) + 
    geom_text(data = tops |> mutate(min=min, max=max, long=min,max),
      aes(y=lat+0.1, x=long+0.1,label='a) Coral cover'), vjust=0, hjust=0) +
    geom_point(data = lgnd.dat, aes(y=y,x=x+0.4, fill=Cat), size=3, shape=21) +
    geom_text(data = lgnd.dat, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0) +
    geom_point(data = towns12 |>
                 filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')),
      aes(y = lat, x = long)) +
    geom_text(data = towns12 |>
                filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')),
      aes(y=lat, x=long, label=town), hjust=1.1, size=3) +
    annotate(geom='text',y=-12,x=145,label=SurveyDates[1], angle=-90, hjust=0.5,vjust=-0.5, size=3) +
    annotate(geom='text',y=-19,x=149.4,label=SurveyDates[2], angle=-31, hjust=0.5,vjust=-0.5, size=3) +
    annotate(geom='text',y=-22.8,x=153.4,label=SurveyDates[3], angle=-73, hjust=0.5,vjust=-0.5, size=3) +
    ## Northern, Central, Southern labels
    annotate(geom='text',y=-12,x=144.5,label=Zones[1], angle=-90, hjust=0.5,vjust=-0.5) +
    annotate(geom='text',y=-18.2,x=147.4,label=Zones[2], angle=-31, hjust=0.5,vjust=-0.5) +
    annotate(geom='text',y=-23,x=153,label=Zones[3], angle=-73, hjust=0.5,vjust=-0.5) +
    scale_x_continuous(breaks=seq(142,154,by=2), position = 'bottom', expand=c(0,0)) +
    scale_y_continuous(breaks=seq(-25,-10,by=2)) +
    coord_sf(ylim=st_bbox(management_sf)[c(2,4)],
       xlim=c(st_bbox(management_sf)[1],
         st_bbox(management_sf)[3]+extra_width)) +
    theme_minimal() +
    theme(legend.position='none', #legend.position=c(0.05,0.10), legend.justification=c(0,0),
          legend.background = element_rect(color='black', fill='white'),
          axis.title=element_blank(),
          panel.border = element_rect(color='black',fill=NA),
          axis.text=element_blank(),
                                        #axis.ticks.length = unit(c(-0.5,-0.9),'lines'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  g.base
}

## ---- Change plot
make_change_plot <- function(manta.mod, management_sf, tops) {
  manta_mod_reefs <- manta.mod |>
    ungroup() |> 
    filter(REPORT_YEAR == final_year) |>
    distinct() |>
    dplyr::select(REEF_NAME)
  lgnd.dat = data.frame(
    x = c(145.2, 145.2),
    y = c(-12, -12.7)
  ) # ff(strt,coefs,length=2)
  lgnd.dat$lab = c("Increase","Decrease")
  lgnd.dat$Cat = factor(c(FALSE,TRUE))
  ## lgnd.dat.size = ff(lgnd.dat[2,1:2],coefs,length=5)
  lgnd.dat.size = ff(c(coefs[1] + ((strt[2]+0.7)*coefs[2]), (strt[2]+0.7)),coefs,length=5)
  lgnd.dat.size$lab = format(seq(2.5,12.5,by=2.5),digits=2)
  lgnd.dat.size$size = seq(2.5,12.5,by=2.5)

  ly.mod <- manta.mod |>
    group_by(REEF_NAME) |>
    summarize(MaxYr=max(REPORT_YEAR, na.rm=TRUE),
      MaxYrCover=Cover[REPORT_YEAR==MaxYr],
      MinYr=max(REPORT_YEAR[REPORT_YEAR<MaxYr], na.rm=TRUE),
      MinYrCover=Cover[REPORT_YEAR==MinYr],
      Longitude=mean(Longitude, na.rm=TRUE),
      Latitude=mean(Latitude, na.rm=TRUE)) |>
    filter(MinYr>=(final_year-2)) |>
    mutate(DiffYr=MaxYr-MinYr,
      Diff = (MaxYrCover - MinYrCover)/DiffYr,
      D=Diff<0) |>
    right_join(manta_mod_reefs) |> filter(!is.na(MaxYr))
  ly.mod |> write.csv(file='../data/processed/FinalYear_coral_change_2023.csv')
    
  rangeLat <- ly.mod |>
    ungroup() |> 
    summarise(across(Latitude, list(Min = min, Max = max),
                     .names = "{.fn}")) |>
    mutate(Diff = Max - Min)
  rangeLong <- ly.mod |>
    ungroup() |> 
    summarise(across(Longitude, list(Min = min, Max = max),
                     .names = "{.fn}")) |>
    mutate(Diff = Max - Min)
  ly1 <- ly.mod |>
    ungroup() |> 
    mutate(Latitude = scales::rescale(Latitude,
                                      to = c(rangeLat$Max, rangeLat$Max + rangeLat$Diff/2)),
           Longitude = scales::rescale(Longitude,
                                       to = c(rangeLong$Max, rangeLong$Max + rangeLong$Diff/2)))

  g.change <-
    ly.mod |>
    ggplot(aes(y = Latitude, x = Longitude)) +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.2, inherit.aes = FALSE) +
    geom_point(aes(y = Latitude,x = Longitude,fill = D, size = abs(Diff*100)),
      alpha = 1,shape = 21,color = 'black', show.legend  =  TRUE)  +
    scale_fill_manual('Change', breaks = c(FALSE,TRUE), labels = c('Increase','Decrease'),
      values = c('green','red'),limits = c(FALSE,TRUE)) +
    geom_text(data = tops |> mutate(min = min, max = max, long = min,max),
      aes(y = lat+0.1, x = long+0.1,label = 'b) Coral change'), vjust = 0, hjust = 0) +
    scale_x_continuous(expand = c(0,0)) +
    geom_point(data = lgnd.dat, aes(y = y,x = x+0.3, fill = Cat),
      shape = 21, size = 3) +
    geom_text(data = lgnd.dat, aes(y = y,x = x+0.3+0.3, label = lab),
      size = 3,hjust = 0, parse = TRUE) +
    geom_point(data = lgnd.dat, aes(y = y,x = x+0.3, fill = Cat),
      shape = 21, size = 3) +
    geom_point(data = lgnd.dat.size, aes(y = y,x = x+0.4, size = size),
      shape = 21) +
    geom_text(data = lgnd.dat.size, aes(y = y,x = x+0.4+0.4, label = lab),
      size = 3,hjust = 0, parse = TRUE) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())    +
    scale_radius()
  g.change
}
add_change_to_plot <- function(g, g.change, management_sf, shiftacross) {
  g.change <- ggplotGrob(g.change)
  g.change.grob.panel = g.change[[1]][[6]]
  m.bb <- sf::st_bbox(management_sf)
  bb.ratio <- diff(m.bb[c(1,3)])/diff(m.bb[c(2,4)])
  ## p.scale.x <- 6
  ## p.scale.y <- p.scale.x/bb.ratio 

  g <- g + annotation_custom(grob = g.change.grob.panel,
                                 xmin = m.bb[1]+shiftacross,
                                 xmax = m.bb[3]+shiftacross,# + p.scale.x,
                                 ymin = -Inf, #- p.scale.y,
                                 ymax = Inf) #bb[2,2] + 0.4) 
  g
}

## ---- COTS
process_cots <- function(cotsAndBleaching) {
  ## Mike supplied the following conversions
  ## COTS = 0: No COTS
  ## COTS>0 & COTS<0.1: No Outbreak
  ## COTS>=0.1 & COTS<0.22: Outbreak watch
  ## COTS>=0.22 & COTS<=1: Incipient Outbreak
  ## COTS>1: Active Outbreak
  ## In addition, if STATUS==RE: Recovering
  cots <- cotsAndBleaching |>
    mutate(OUTBREAK.CAT = case_when(
      AvgOfMEAN_COTS == 0 ~ "No COTS",
      AvgOfMEAN_COTS > 0 & AvgOfMEAN_COTS < 0.1 ~ "No Outbreak",
      AvgOfMEAN_COTS >= 0.1 & AvgOfMEAN_COTS < 0.22 ~ "Outbreak Watch",
      AvgOfMEAN_COTS >= 0.22 & AvgOfMEAN_COTS <= 1 ~ "Established Outbreak", #' Incipient Outbreak',
      AvgOfMEAN_COTS > 1 ~ "Severe Outbreak", #' Active Outbreak',
      STATUS == "RE" ~ "Recovering"
    )) |>
    mutate(
      OUTBREAK.CAT = ifelse(OUTBREAK.CAT == "Outbreak watch",
        "Outbreak Watch", OUTBREAK.CAT),
      OUTBREAK.CAT = factor(OUTBREAK.CAT,
        levels = c("No COTS", "No Outbreak", "Outbreak Watch",
          "Recovering", "Established Outbreak", "Severe Outbreak")
      )
    )
  ## Actually, it turns out that 'Outbreak Watch' should be 'Potential Outbreak'
  cots <- cots |>
    mutate(
      OUTBREAK.CAT = as.character(OUTBREAK.CAT),
      OUTBREAK.CAT = ifelse(OUTBREAK.CAT == "Outbreak Watch",
        "Potential Outbreak", OUTBREAK.CAT
      ),
      OUTBREAK.CAT = factor(OUTBREAK.CAT,
        levels = c("No COTS", "No Outbreak", "Potential Outbreak",
          "Recovering", "Established Outbreak", "Severe Outbreak")
      )
    )
  cots
}
make_cots_plot <- function(cots, coefs, tops) {
  lgnd.dat <- ff(c(coefs[1] + (strt[2]+0.7)*coefs[2], strt[2]+0.7),coefs,length=6)
  lgnd.dat$lab <- levels(cots$OUTBREAK.CAT)
  lgnd.dat$Cat <- factor(lgnd.dat$lab)
  g.cots = cots |>
    dplyr::rename(Latitude = REEF_LAT, Longitude = REEF_LONG) |>
    ggplot(aes(y = Latitude, x = Longitude)) +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.2, inherit.aes = FALSE) +
    geom_point(aes(y = Latitude, x = Longitude, fill = OUTBREAK.CAT),
      size = 3, alpha = 1, shape = 21, color = "black", show.legend = TRUE
    ) +
    scale_fill_manual("COTS",
      breaks = levels(cots$OUTBREAK.CAT),
      labels = levels(cots$OUTBREAK.CAT),
      values = rev(trafficLightColors(6)),
      limits = levels(cots$OUTBREAK.CAT)
    ) +
    geom_text(
      data = tops |>
        mutate(min = min, max = max, long = min, max),
      aes(y = lat + 0.1, x = long + 0.1, label = "c) COTS"),
      vjust = 0, hjust = 0
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_point(data = lgnd.dat, aes(y = y, x = x + 0.5, fill = Cat), shape = 21, size = 3) +
    geom_text(
      data = lgnd.dat, aes(y = y, x = x + 0.5 + 0.3, label = lab),
      size = 3, hjust = 0, parse = FALSE
    ) +
    geom_point(
      data = lgnd.dat, aes(y = y, x = x + 0.5, fill = Cat),
      shape = 21, size = 3
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  g.cots
  
}

## ---- Bleaching
process_bleaching <- function(cotsAndBleaching, cots) {
  bleaching <- cotsAndBleaching |>
    mutate(`Bleaching Aes` = Ave_perc_bleach_cat)                  # 2023
  bleaching <- bleaching |>
    mutate(
      CAT = ifelse(is.na(`Bleaching Aes`), 0,
        ifelse(`Bleaching Aes` == 0, 0,
          ifelse(`Bleaching Aes` %in% c("0+", "1-", "1", "1+"), "<10%",
            ifelse(`Bleaching Aes` %in% c("2-", "2", "2+"), "10-30%",
              ifelse(`Bleaching Aes` %in% c("3-", "3", "3+"), "30-60%", ">60%")
            )
          )
        )
      ),
      CAT = factor(CAT, levels = c("0", "<10%", "10-30%", "30-60%", ">60%"))
    ) |>
    left_join(cots |>
                dplyr::select(REEF_LAT, REEF_LONG, REEF_NAME) |>
                distinct())
  bleaching
}
make_bleaching_plot <- function(bleaching, coefs, tops, management_sf) {
  bleaching <- bleaching |>
    mutate(Date = as.Date(`in-water_sample_date`, format = "%d/%m/%Y"),
           cDate = ifelse(Date > as.Date(paste0(final_year, "-01-01")), "After", "Before"),
           cDate = factor(cDate, levels = c('Before', 'After')) )
  levels(bleaching$CAT) <-  c('0%', '>0% - 10%', '>10% - 30%', '>30% - 60%', '>60%') 
  lgnd.dat <- ff(c(coefs[1] + (strt[2]+0.7)*coefs[2], (strt[2]+0.7)), coefs, length = 5)
  lgnd.dat$lab <- levels(bleaching$CAT)
  lgnd.dat$Cat <- factor(lgnd.dat$lab)
  g.bleaching <-
    bleaching |>
    dplyr::rename(Latitude = REEF_LAT, Longitude = REEF_LONG) |>
    ggplot(aes(y = Latitude, x = Longitude)) +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.2, inherit.aes = FALSE) +
    geom_point(aes(
      y = Latitude, x = Longitude, fill = CAT,
      shape = cDate
    ), size = 3, alpha = 1, color = "black", show.legend = TRUE) +
    scale_fill_manual("Bleaching",
      breaks = c(levels(bleaching$CAT)),
      labels = c(levels(bleaching$CAT)),
      values = c(rev(trafficLightColors(5))),
      limits = c(levels(bleaching$CAT))
    ) +
    scale_shape_manual("Survey Date",
      breaks = c("Before", "After"),
      labels = c(paste0("Prior to 1st Jan ", final_year), paste0("Post 1st Jan ", final_year)),
      values = c(24, 21)
    ) +
    geom_text(
      data = tops |>
        mutate(min = min, max = max, long = min, max),
      aes(y = lat + 0.1, x = long + 0.1, label = "d) Bleaching - in-water surveys"),
      vjust = 0, hjust = 0
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_point(
      data = lgnd.dat, aes(y = y, x = x + 0.5, fill = Cat),
      shape = 21, size = 3
    ) +
    geom_text(
      data = lgnd.dat, aes(y = y, x = x + 0.5 + 0.3, label = lab),
      size = 3, hjust = 0, parse = FALSE
    ) +
    geom_point(
      data = lgnd.dat, aes(y = y, x = x + 0.5, fill = Cat),
      shape = 21, size = 3
    ) +
    annotate(geom = "point", y = c(-11.5), x = c(145.3), shape = 24, size = 3) +
    annotate(
      geom = "text", y = c(-11.5), x = c(145.6),
      label = paste0("Prior to 1st Jan ", final_year), hjust = 0, size = 3
    ) +
    annotate(
      geom = "point", y = c(-12.2), x = c(145.3),
      shape = 21, size = 3
    ) +
    annotate(
      geom = "text", y = c(-12.2), x = c(145.6),
      label = paste0("Post 1st Jan ", final_year), hjust = 0, size = 3
    ) +
    ## geom_point(data = NULL, aes(y = c(-12, -13), x = c(145, 145), shape = c(1, 2))) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  g.bleaching
}

## ---- Aerial Bleaching
process_aerial_bleaching <- function(cotsAndBleaching) {
  bleaching.aerial <- cotsAndBleaching |>
     mutate(`Bleaching Aerial` = Aerial_bleaching)
  bleaching.aerial <- bleaching.aerial |>
     mutate(
       CAT = ifelse(is.na(`Bleaching Aerial`), NA,
         ifelse(`Bleaching Aerial` == 0, 0,
           ifelse(`Bleaching Aerial` %in% c("0+", "1-", "1", "1+"), "<10%",
             ifelse(`Bleaching Aerial` %in% c("2-", "2", "2+"), "10-30%",
               ifelse(`Bleaching Aerial` %in% c("3-", "3", "3+"), "30-60%", ">60%")
             )
           )
         )
       ),
       ## CAT = factor(CAT, levels=c(NA,'0','<10%','10-30%','30-50%','>50%'))) |>
       CAT = factor(CAT, levels = c("0", "<10%", "10-30%", "30-60%", ">60%"))
     ) |>
     left_join(cots |>
     dplyr::select(REEF_LAT, REEF_LONG, REEF_NAME) |>
     distinct())
  bleaching.aerial1 <- bleaching.aerial |>
     filter(!is.na(CAT)) |>
     droplevels()
  bleaching.aerial1
}
make_nice_date_range <- function(bleaching.aerial_date_range) {
  bleaching.aerial_date_range <- bleaching.aerial |>
    group_by(Region) |>
    summarise(min = min(Aerial_date, na.rm = TRUE), max = max(Aerial_date, na.rm = TRUE))
  bleaching.aerial_date_range |>
    mutate(survey_date = ifelse(month(min) == month(max),
      paste0(day(min), "-", day(max), " ", month(min, label = TRUE), " ", final_year),
      paste0(
        day(min), " ", month(min, label = TRUE), " ", final_year, " - ",
        day(max), " ", month(max, label = TRUE), " ", final_year
      )
    )) |>
    dplyr::select(Region, survey_date)
}
make_aerial_bleaching_plot <- function(bleaching.aerial, coefs, tops, management_sf) {
  survey_dates <- make_nice_date_range(bleacing.aerial) |>
    mutate(Region = factor(Region, levels = c("Northern", "Central", "Southern"))) |>
    arrange(Region)
  SurveyDates <- paste0("Suvery dates\n(", survey_dates |> pull(survey_date), ")")

  levels(bleaching.aerial$CAT) <-  c('0','>0% - 10%', '>10% - 30%', '>30% - 60%', '>60%') 

  ## lgnd.dat <- ff(c(coefs[1] + (strt[2]*coefs[2]), strt[2]),coefs,length=5)
  lgnd.dat = ff(c(coefs[1] + ((strt[2]+0.7)*coefs[2]), (strt[2]+0.7)),coefs,length=5)
  lgnd.dat$lab = levels(bleaching.aerial$CAT)
  lgnd.dat$Cat=factor(lgnd.dat$lab)
  g.bleaching.aerial <-
    bleaching.aerial |>
    dplyr::rename(Latitude = REEF_LAT, Longitude = REEF_LONG) |>
    ggplot(aes(y = Latitude, x = Longitude)) +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.2, inherit.aes = FALSE) +
    annotate(
      geom = "text", y = -12.8, x = 145.2, label = SurveyDates[1],
      angle = -0, hjust = 0, vjust = -0.5, size = 3
    ) +
    annotate(
      geom = "text", y = -19, x = 149.4, label = SurveyDates[2],
      angle = -31, hjust = 0.5, vjust = -0.5, size = 3
    ) +
    annotate(
      geom = "text", y = -22.8, x = 153.4, label = SurveyDates[3],
      angle = -73, hjust = 0.5, vjust = -0.5, size = 3
    ) +
    geom_point(aes(y = Latitude, x = Longitude, fill = CAT),
      size = 3, alpha = 1, shape = 21, color = "black", show.legend = TRUE) +
    scale_fill_manual("Bleaching",
      breaks = c(levels(bleaching.aerial$CAT)),
      labels = c(levels(bleaching.aerial$CAT)),
      values = c(rev(trafficLightColors(5))),
      limits = c(levels(bleaching.aerial$CAT))
    ) +
    geom_text(data = tops |>
                mutate(min = min, max = max, long = min, max),
      aes(y = lat + 0.1, x = long, label = "e) Bleaching - aerial surveys"),
      vjust = 0, hjust = 0) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_point(data = lgnd.dat, aes(y = y + 0.4, x = x + 0.25,
      fill = Cat), shape = 21, size = 3) +
    geom_text(data = lgnd.dat, aes(y = y + 0.4, x = x + 0.25 + 0.3,
      label = lab), size = 3, hjust = 0, parse = FALSE) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  g.bleaching.aerial 
}

## ---- Oz inset
make_oz_inset <- function(qld_sf) {
  aus <- rnaturalearth::ne_countries(scale = "medium", country = "Australia") |>
    st_crop(c(xmin = 110, ymin = -45, xmax = 155, ymax = -5))

  qld_sf <-
    qld_sf |>
    st_make_valid() |> 
    st_crop(
      c(
        st_bbox(management_sf)[1],
        st_bbox(management_sf)[2] - (st_bbox(management_sf)[4] - st_bbox(management_sf)[2]) * 0.1,
        st_bbox(management_sf)[3],
        st_bbox(management_sf)[4]
      ))
  aus <-
    aus |>
    ggplot() +
    geom_sf(fill = "white", colour = "black") +
    geom_sf(data = qld_sf, fill = "grey", colour = "black") +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.2, inherit.aes = FALSE) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          panel.background = element_rect(color='black'))
  aus 
}
add_oz_to_plot <- function(g, g.oz, management_sf, extra_width, inset_width = 5, xmin, ymin) {
  aus <- rnaturalearth::ne_countries(scale = "medium", country = "Australia") |>
    st_crop(c(xmin = 110, ymin = -45, xmax = 155, ymax = -5))
  bb.ratio <- (st_bbox(aus)[3] - st_bbox(aus)[1]) / (st_bbox(aus)[4] - st_bbox(aus)[2])

  bb.xmin <- st_bbox(management_sf)[1]
  bb.xmax <- st_bbox(management_sf)[3] + extra_width
  bb.ymin <- st_bbox(management_sf)[2]
  bb.ymax <- st_bbox(management_sf)[4]
  
  g.oz.grob=ggplotGrob(g.oz)
  g = g + annotation_custom(
    grob = g.oz.grob, xmax = Inf, ymax = Inf,
    xmin = bb.xmax - inset_width,
    ymin = bb.ymax - (inset_width/bb.ratio)
  )
  g
}

## ---- Scalebar and N arrow
add_scale_and_arrow_to_plot <- function(g) {
  gc.1=c(143.5,-24)
  gc.2=gcDestination(lon=gc.1[1], lat=gc.1[2], bearing=90, dist=500, model='WGS84', Vincenty=FALSE)
  g = g + ggsn:::scalebar(x.min=gc.1[1],y.min=gc.1[2],x.max=gc.2[1],y.max=gc.1[2]+0.5,
    dist=250, model='WGS84',st.size=3, height=0.5,st.dist=0.5, dist_unit = 'km',
    transform=TRUE) +
    ggsn:::north(x.min=mean(c(gc.1[1],gc.2[1]))-0.5, x.max=mean(c(gc.1[1],gc.2[1]))+0.5, y.min=gc.1[2]+0.5,y.max=gc.1[2]+2, scale=1) 
  g1=g
  g1
}


gcDestination <- function(lon, lat, bearing, dist, dist.units = "km",
    model=NULL, Vincenty=FALSE) {
 # lat, lon : lattitude and longitude in decimal degrees
 # bearing : bearing from 0 to 360 degrees
 # dist : distance travelled
 # dist.units : units of distance "km" (kilometers), "nm" (nautical 
 # miles), "mi" (statute miles)
 # model : choice of ellipsoid model ("WGS84", "GRS80", "Airy", 
 # "International", "Clarke", "GRS67")

    if (!is.numeric(lon)) stop("lon not numeric")
    if (!is.numeric(lat)) stop("lat not numeric")
    if (!is.numeric(bearing)) stop("bearing not numeric")
    if (!is.numeric(dist)) stop("dist not numeric")

    if (length(lon) != length(lat)) stop("lon and lat differ in length")
    if (length(bearing) > 1L && length(lon) > 1L) stop("length mismatch")
    if (length(bearing) > 1L && length(dist) > 1L) stop("length mismatch")

    as.radians <- function(degrees) degrees * pi / 180
    as.degrees <- function(radians) radians * 180 / pi
    as.bearing <- function(radians) (as.degrees(radians) + 360) %% 360

    ellipsoid <- function(model = "WGS84") {
        switch(model,
        WGS84 = c(a = 6378137, b = 6356752.3142, f = 1 / 298.257223563),
        GRS80 = c(a = 6378137, b = 6356752.3141, f = 1 / 298.257222101),
        Airy = c(a = 6377563.396, b = 6356256.909, f = 1 / 299.3249646),
        International = c(a = 6378888, b = 6356911.946, f = 1 / 297),
        Clarke = c(a = 6378249.145, b = 6356514.86955, f = 1 / 293.465),
        GRS67 = c(a = 6378160, b = 6356774.719, f = 1 / 298.25),
        c(a = NA, b = NA, f = NA)
    )}
    
    dist <- switch(dist.units,
        km = dist,
        nm = dist * 1.852,
        mi = dist * 1.609344
    )
    lat <- as.radians(lat)
    lon <- as.radians(lon)
    bearing <- as.radians(bearing)

    if (is.null(model)) {
 # Code adapted from JavaScript by Chris Veness 
 # (scripts@movable-type.co.uk) at
 # http://www.movable-type.co.uk/scripts/latlong.html#ellipsoid
 #   originally from Ed Williams' Aviation Formulary, 
 # http://williams.best.vwh.net/avform.htm
        radius <- 6371
        psi <- dist / radius
        lat2 <- asin(sin(lat) * cos(psi) +  cos(lat) * sin(psi) * cos(bearing))
        lon2 <- lon + atan2(sin(bearing) * sin(psi) * cos(lat), cos(psi) - 
            sin(lat) * sin(lat2))
        if (any(is.nan(lat2)) || any(is.nan(lon2))) warning("Out of range values")
        return(cbind(long=as.degrees(lon2), lat=as.degrees(lat2)))
    }

    ellips <- ellipsoid(model)
    if (is.na(ellips["a"])) stop("no such ellipsoid model")
    if (Vincenty) {
 # Code adapted from JavaScript by Chris Veness 
 # (scripts@movable-type.co.uk) at
 # http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
 # Original reference (http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf):
 #   Vincenty, T. 1975.  Direct and inverse solutions of geodesics on 
 # the ellipsoid with application of nested equations.
 #      Survey Review 22(176):88-93
        dist <- dist * 1000
        sin.alpha1 <- sin(bearing)
        cos.alpha1 <- cos(bearing)
        tan.u1 <- (1 - ellips["f"]) * tan(lat)
        cos.u1 <- 1 / sqrt(1 + (tan.u1 ^ 2))
        sin.u1 <- tan.u1 * cos.u1
        sigma1 <- atan2(tan.u1, cos.alpha1)
        sin.alpha <- cos.u1 * sin.alpha1
        cos.sq.alpha <- 1 - (sin.alpha ^ 2)
        u.sq <- cos.sq.alpha * ((ellips["a"] ^ 2) - (ellips["b"] ^ 2)) /
            (ellips["b"] ^ 2)
        cap.A <- 1 + u.sq / 16384 * (4096 + u.sq * (-768 + u.sq * (320 - 
            175 * u.sq)))
        cap.B <- u.sq / 1024 * (256 + u.sq * (-128 + u.sq * (74 - 47 * u.sq)))

        sigma <- dist / (ellips["b"] * cap.A)
        sigma.p <- 2 * pi
        cos.2.sigma.m <- cos(2 * sigma1 + sigma)
        while(any(abs(sigma - sigma.p) > 1e-12)) {
            cos.2.sigma.m <- cos(2 * sigma1 + sigma)
            sin.sigma <- sin(sigma)
            cos.sigma <- cos(sigma)
            delta.sigma <- cap.B * sin.sigma * (cos.2.sigma.m + cap.B / 4 * 
                (cos.sigma *
                (-1 + 2 * cos.2.sigma.m ^ 2) - cap.B / 6 * cos.2.sigma.m *
                (-3 + 4 * sin.sigma ^ 2) * (-3 + 4 * cos.2.sigma.m ^ 2)))
            sigma.p <- sigma
            sigma <- dist / (ellips["a"] * cap.A) + delta.sigma
        }
        tmp <- sin.u1 * sin.sigma - cos.u1 * cos.sigma * cos.alpha1
        lat2 <- atan2(sin.u1 * cos.sigma + cos.u1 * sin.sigma * cos.alpha1,
            (1 - ellips["f"]) * sqrt(sin.alpha ^ 2 + tmp ^ 2))
        lambda <- atan2(sin.sigma * sin.alpha1, cos.u1 * cos.sigma - sin.u1 * 
            sin.sigma * cos.alpha1)
        cap.C <- ellips["f"] / 16 * cos.sq.alpha * (4 + ellips["f"] * 
            (ellips["f"] - 3 * cos.sq.alpha))
        cap.L <- lambda - (1 - cap.C) * ellips["f"] * sin.alpha *
            (sigma + cap.C * sin.sigma * (cos.2.sigma.m + cap.C * cos.sigma * 
            (-1 + 2 * cos.2.sigma.m ^ 2)))
        lat2 <- as.degrees(lat2)
        lon2 <- as.degrees(lon + cap.L)
    } else {
 # Code adapted from JavaScript by Larry Bogan (larry@go.ednet.ns.ca) 
 # at http://www.go.ednet.ns.ca/~larry/bsc/jslatlng.html
        e <- 0.08181922
        radius <- (ellips["a"] / 1000) * (1 - e^2) / ((1 - e^2 * 
            sin(lat)^2)^1.5)
        psi <- dist / radius
        phi <- pi / 2 - lat
        arc.cos <- cos(psi) * cos(phi) + sin(psi) * sin(phi) * cos(bearing)
        lat2 <- as.degrees((pi / 2) - acos(arc.cos))
        arc.sin <- sin(bearing) * sin(psi) / sin(phi)
        lon2 <- as.degrees(lon + asin(arc.sin))
    }
    return(cbind(long=lon2, lat=lat2))
}

## ---- Watermarks
add_watermark_to_plot <- function(g, water_mark, management_sf, extra_width, oz_width, inset_width) {
  aus <- rnaturalearth::ne_countries(scale = "medium", country = "Australia") |>
    st_crop(c(xmin = 110, ymin = -45, xmax = 155, ymax = -5))
  bb.ratio <- (st_bbox(aus)[3] - st_bbox(aus)[1]) / (st_bbox(aus)[4] - st_bbox(aus)[2])

  bb.xmin <- st_bbox(management_sf)[1]
  bb.xmax <- st_bbox(management_sf)[3] + extra_width
  bb.ymin <- st_bbox(management_sf)[2]
  bb.ymax <- st_bbox(management_sf)[4]
  im.width <- image_info(water_mark)$width
  im.height <- image_info(water_mark)$height
  im.ratio <- im.width / im.height

  im.width <- (inset_width / (bb.xmax - bb.xmin))
  im.x <- (1 - (oz_width / (bb.xmax - bb.xmin))/2)
  ## im.x <- 1 - im.width
  im.y <- (1 - im.width) * im.ratio
  
  g + 
    annotation_custom(rasterGrob(water_mark,
                                 ## x=unit(0.90, 'npc'),
                                 x=unit(im.x, 'npc'),
                                 ## x=unit(0.89, 'npc'),
                                 ## y=unit(0.65, 'npc'),
                                 y=unit(im.y, 'npc'),
                                 vjust = 1,
                                 hjust = 0.5,
                                 width=unit(im.width, 'npc')))
                                 ## width=unit(0.1, 'npc')))
}






