#######################################################################
## The following function queries the                                ##
## /export/project/monitoring-dashboard/data/dashboard.sqlite db     ##
## to assess what reefs are represented by data and where those data ##
## are stored.                                                       ##
## Returns
##   a data.frame containing the raw data for the modelled reefs     ##
#######################################################################
CoralTrends_get_raw_data_via_db_and_file <- function() {
  data <- get_candidates(tab_name = "not_needed", data_type = "manta", scale = "reef") |>
    ## filter(selected_flag == 1) |>
    filter(selected_flag == TRUE | selected_flag == 1) |>
    mutate(nm = paste0(config_$model_path,
      paste(data_type, data_scale, domain_name,
        group,
        family_type,
        reef_zone, depth,
        shelf,
        ## selectors_list$response_selector,
        "Cover",
        sub_model,
        "raw_data.rds", sep = "_"))) |>
    mutate(raw = map(.x = nm,
      .f = ~ {
        if (file.exists(.x)) {
          readRDS(.x)
        } else NULL
      })) |>
    dplyr::select(nm, raw) |>
    unnest("raw")
  return(data)
}

#######################################################################
## The following function removes "_NA" strings from paths           ##
## Returns                                                           ##
##   a data.frame containing the raw data for the modelled reefs     ##
#######################################################################
CoralTrends_tidy_raw_data <- function(raw_manta) {
  raw_manta <- raw_manta |>
    mutate(REEF = str_replace(REEF, " NA", "")) |>
    group_by(REEF) |>
    mutate(
      LATITUDE = mean(LATITUDE),
      LONGITUDE = mean(LONGITUDE)
    ) |>
  filter(SECTOR != "TS")
  return(raw_manta)
}

#######################################################################
## The following function takes a vector of paths and returns a      ##
## single path that represents the common elements across all paths  ##
## Returns                                                           ##
##   a character string path                                         ##
#######################################################################
find_common_pattern <- function(strings) {
  if (length(strings) == 0) return("")
  # Split each string into parts based on "_"
  split_strings <- strsplit(strings, "_", fixed = TRUE)

  # Find the length of the shortest split (to avoid index errors)
  min_length <- min(sapply(split_strings, length))

  # Compare elements at each position
  result <- sapply(seq_len(min_length), function(i) {
    elements <- sapply(split_strings, `[`, i)
    if (all(elements == elements[1])) {
      return(elements[1])  # Keep common elements
    } else {
      return(" ")  # Replace differing elements with "*"
    }
  })

  # Collapse into a single string
  paste(result, collapse = "_")
}

#######################################################################
## The following function queries the                                ##
## /export/project/monitoring-dashboard/data/dashboard.sqlite db     ##
## to assess what reefs have annual estimates or annual comparisons  ##
## summarised and where those summaries are stored                   ##
## Returns                                                           ##
##   a list comprising:                                              ##
##      - data.frame containing a path to the stored data            ##
##      - a path with a common name (e.g. reef name removed)         ##
#######################################################################
CoralTrends_get_annual_data_via_db_and_file <- function(type = "year_sum") {
  data <- get_candidates(tab_name = "not_needed", data_type = "manta", scale = "reef") |>
    filter(selected_flag == TRUE | selected_flag == 1) |>
    mutate(nm = paste0(config_$model_path,
      paste(data_type, data_scale, domain_name,
        group,
        family_type, reef_zone, depth,
        shelf, model_type, sub_model,
        paste0(type, ".rds"), sep = "_"))) |>
    dplyr::select(data_type, data_scale, domain_name,
      group, family_type, reef_zone, depth,
      shelf, model_type, sub_model, nm)
  nm <- find_common_pattern(data$nm)
  return(list(data = data, nm = nm))
}

#######################################################################
## The following function reads each of the annual summaries into a  ##
## single data frame                                                 ##
## Returns                                                           ##
##   a data frame of all reefs by thier measured years               ##
#######################################################################
CoralTrends_tidy_annual_manta_data <- function(annual_data,
                                               raw_data,
                                               final_year = 2000,
                                               type = "year_sum") {
  annual_manta <- annual_data$data |>
    mutate(dat =  map2(.x = nm, .y = reef_zone,
      .f = ~ {
        if (!file.exists(.x)) return(NULL)
        x <- readRDS(.x)
        return(x)
      })) |>
    unnest(dat) |>
    dplyr::select(-shelf)
  if (type == "year_sum") {
    annual_manta <- annual_manta |>
      ## left_join(raw_manta |>
      inner_join(raw_manta |>
                  dplyr::select(domain_name = REEF,
                    Latitude = LATITUDE,
                    Longitude = LONGITUDE) |>
                  distinct(),
        by = "domain_name")
  } else {
    annual_manta <- annual_manta |>
      ## left_join(
      inner_join(
        raw_manta |>
          dplyr::select(domain_name = REEF, Latitude = LATITUDE, Longitude = LONGITUDE) |>
          distinct(),
        by = "domain_name"
      ) |>
      separate(YearComp, into = c("Yr2", "Yr1"), sep = "-", remove = FALSE) |>
      filter(Yr2 == final_year) %>%
      mutate(
        Yr2 = as.numeric(Yr2),
        Yr1 = as.numeric(Yr1)
      ) %>%
      arrange(domain_name, desc(YearComp)) %>%
      group_by(domain_name) %>%
      slice(1) %>%
      filter(Yr1 >= (final_year - 2)) %>%
      mutate(
        DiffYr = Yr2 - Yr1,
        Diff = mean / DiffYr,
        D = ifelse(Pl > 0.95, "decrease", ifelse(Pg > 0.95, "increase", "neutral")),
        D = factor(D, levels = c("increase", "neutral", "decrease"))
      ) |>
      mutate(REPORT_YEAR = Yr2)

  }
  return(annual_manta)
}

#######################################################################
## The following function combine the coral cover and coral change   ##
## data frames into a single data frame                              ##
## Returns                                                           ##
##   a data frame that combines coral cover and coral change         ##
#######################################################################
CoralTrends_combine_cover_and_change <- function(annual_manta,
                                                 annual_yearcomp_manta,
                                                 final_year = 2000) {
  coral_cover_and_change <-
    annual_manta |>
    filter(REPORT_YEAR == final_year) %>%
    mutate(variable = "Cover") |>
    dplyr::select(REPORT_YEAR,
      REEF_NAME = domain_name, variable,
      Mean = mean, Median = median, lower, upper, Latitude, Longitude
    ) |>
    bind_rows(
      annual_yearcomp_manta |>
        filter(variable == "value") |>
        mutate(variable = "AnnualDiff") |>
        dplyr::select(REPORT_YEAR,
          REEF_NAME = domain_name, variable,
          Mean = Diff, Median = Diff, lower, upper, Latitude, Longitude, D
        )
    ) |>
    mutate(latitude = Latitude, longitude = Longitude) |>
    arrange(REPORT_YEAR, REEF_NAME, variable)
  return(coral_cover_and_change)
}

CoralTrends_process_manta_mod_data__ <- function(coral_cover_and_change) {
  coral_cover_and_change |>
    mutate(Latitude = latitude, Longitude = longitude, Cover = Median)
}

#
# Convert GBR management_sf into a data.frame
get_gbr_df <- function(management_sf) {
  management_sf |>
    st_coordinates() |>
    as_tibble() |>
    rename(long = X, lat = Y)
}

## Get the lat/long coordinates for the top of the GBR polygons
get_gbr_tops <- function(management_df) {
  ## Postion of labels above GBR polygons
  management_df |>
    filter(lat > max(lat) - 0.1) |>
    summarize(min = min(long), max = max(long), lat = max(lat))
}

## Get the slope of the GBR
get_management_coefs <- function(management_df) {
  dat <- management_df |>
    filter(lat < -15 & lat > -16, long > 145.5)
  coefs <- coef(lm(long ~ lat, data = dat))
  coefs
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

## Define the position of the "Northern", "Central" and "Southern" labels
get_lgnd_pos <- function(coefs) {
  strt = c(145.9,-14.5)

  lgnd.dat <- ff(c(coefs[1] + (strt[2]*coefs[2]), strt[2]),coefs,length=5)
  lgnd.dat$lab <- rev(c('>0% - 10%','>10% - 30%','>30% - 50%','>50% - 75%','>75% - 100%'))
  lgnd.dat$Cat <- factor(5:1)
  lgnd.dat
}

## Define the traffic light colors to use for coral change fill
trafficLightColors <- colorRampPalette(c('#ED1C24','#F47721','#F0C918','#B0D235','#00734D'))

## Define the towns to plot on the Queensland map
towns12 <-
  structure(list(town = c("Cooktown", "Cairns", "Townsville", "Bowen",
    "Proserpine", "Mackay", "Rockhampton", "Gladstone", "Port Douglas",
    "Ayr", "Innisfail", "Cardwell"), long = c(145.229629516602, 145.740203857422,
      146.78092956543, 148.190048217773, 148.696716308594, 149.112060546875,
      150.409729003906, 151.094924926758, 145.405410766602, 147.405349731445,
      146.016479492188, 146.024261474609), lat = c(-15.4718799591065,
        -16.9359798431396, -19.2582454681396, -19.9939575195313, -20.311653137207,
        -21.1522636413574, -23.3642120361328, -23.8742027282715, -16.4558181762695,
        -19.5710220336914, -17.516508102417, -18.2670269012451)), .Names = c("town",
          "long", "lat"), row.names = c("Cooktown", "Cairns", "Townsville",
            "Bowen", "Proserpine", "Mackay", "Rockhampton", "Gladstone",
            "Port Douglas", "Ayr", "Innisfail", "Cardwell"), class = "data.frame", col = 1, pch = 21, bg = "red", cex = 1, txt.text = c("Cooktown",
              "Cairns", "Townsville", "Bowen", "Proserpine", "Mackay", "Rockhampton",
              "Gladstone", "Port Douglas", "Ayr", "Innisfail", "Cardwell"), txt.col = 1, txt.cex = 1, txt.pos = 2)

CoralTrends_get_cots_and_bleaching <- function() {
  ## cat(paste0(sql_path, "cots_and_bleaching.sql"), "\n")
  ## cat(paste0(data_path, "cots_and_bleaching.csv"), "\n")
  ## cat(file.exists(paste0(sql_path, "cots_and_bleaching.sql")))
  system2("java",
    args = c("-jar",
      "/export/project/monitoring-dashboard/dev/dbExport.jar",
      shQuote(paste0(sql_path, "cots_and_bleaching.sql")),
      shQuote(paste0(data_path, "cots_and_bleaching.csv")),
      "reef",
      "reefmon",
      paste0(">> '", log_file, "' 2>&1")),
    wait = TRUE,
    stderr = sub("dashboard.log$", "dashboard_error.log", config_$dashboard_log)
  )
  cots_and_bleaching <- read_csv(paste0(data_path, "cots_and_bleaching.csv"))

  cots_and_bleaching <-
    cots_and_bleaching |>
    mutate(sample_date = as.POSIXct(SAMPLE_DATE, format = "%d-%b-%Y %H:%M:%S")) |>
    mutate(sample_date = format(sample_date, "%d/%m/%Y")) |>
    group_by(AIMS_REEF_NAME, REEF_LAT, REEF_LONG, REPORT_YEAR) |>
    summarise(
      number_of_tows = length(TOW_SEQ_NO),
      AvgOfMEAN_COTS = mean(COT_COUNT, na.rm = TRUE),
      cots_per_tow = AvgOfMEAN_COTS / number_of_tows,
      MaxOfBLEACHING_PERHC = max(BLEACHING_PERHC, na.rm = TRUE),
      `in-water_sample_date` = max(sample_date, na.rm = TRUE),
      .groups = "drop") |>
    mutate(ltmp_bleach_cat = str_replace(MaxOfBLEACHING_PERHC, "U|L", "")) |>
    rename(REEF = AIMS_REEF_NAME) |>
    left_join(
      CoralTrends_get_raw_data_via_db_and_file() |>
        mutate(REEF = str_replace(REEF, " NA", "")) |>
        dplyr::select(Region, REEF) |>
        distinct(),
      by = "REEF"
    ) |>
    (\(x) {
      if (!"aerial_bleach" %in% names(x)) {
        dplyr::mutate(x, aerial_bleach = NA)
      } else {
        x
      }
    })() |>
     (\(x) {
       if (!"aerial_date" %in% names(x)) {
         dplyr::mutate(x, aerial_date = NA)
       } else {
         x
       }
     })() |>
      mutate(STATUS = case_when(
        cots_per_tow == 0 ~ "NC",
        cots_per_tow > 0 & cots_per_tow <= 0.1 ~ "NO",
        cots_per_tow > 0.1 & cots_per_tow <= 0.22 ~ "PO",
        cots_per_tow > 0.22 & cots_per_tow <= 1 ~ "EO",
        cots_per_tow > 1 ~ "SO"
      )) |>
      filter(REPORT_YEAR == final_year)

  cots_and_bleaching <-
    cots_and_bleaching |>
    mutate(Ave_perc_bleach_cat = as.character(ltmp_bleach_cat),
      Aerial_bleaching = as.character(aerial_bleach),
      Aerial_date = ifelse(is.na(aerial_date),
        NA,
        as.Date(aerial_date = format("%d/%m/%Y"))))

  return(cots_and_bleaching)
}

## COTS
process_cots <- function(cots_and_bleaching) {
  ## Mike supplied the following conversions
  ## COTS = 0: No COTS
  ## COTS>0 & COTS<0.1: No Outbreak
  ## COTS>=0.1 & COTS<0.22: Outbreak watch
  ## COTS>=0.22 & COTS<=1: Incipient Outbreak
  ## COTS>1: Active Outbreak
  ## In addition, if STATUS==RE: Recovering
  cots <- cots_and_bleaching |>
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
make_cots_plot <- function(cots, coefs, tops, lgnd.dat, management_sf) {
  strt <- c(145.9,-14.5)
  lgnd.dat <- ff(c(coefs[1] + (strt[2]+0.7)*coefs[2], strt[2]+0.7),coefs,length=6)
  lgnd.dat$lab <- levels(cots$OUTBREAK.CAT)
  lgnd.dat$Cat <- factor(lgnd.dat$lab)
  g.cots <- cots |>
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

## Bleaching
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
                ## dplyr::select(REEF_LAT, REEF_LONG, REEF_NAME) |>
                dplyr::select(REEF_LAT, REEF_LONG, REEF) |>
                distinct())
  bleaching
}

make_bleaching_plot <- function(bleaching, coefs, tops, lnd.dat, management_sf) {
  strt <- c(145.9,-14.5)
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

## Aerial Bleaching
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
       CAT = factor(CAT, levels = c("0%", "<10%", "10-30%", "30-60%", ">60%"))
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
    mutate(survey_date =
             ifelse(is.infinite(min),
               "Not surveyed",
               ifelse(month(min) == month(max),
                 ifelse(day(min) == day(max),
                   paste0(day(min), " ", month(min, label = TRUE), " ", final_year),
                   paste0(day(min), "-", day(max), " ", month(min, label = TRUE), " ", final_year)),
                 paste0(
                   day(min), " ", month(min, label = TRUE), " ", final_year, " - ",
                   day(max), " ", month(max, label = TRUE), " ", final_year
                 )
               )
             )) |>
    dplyr::select(Region, survey_date)
}

make_aerial_bleaching_plot <- function(bleaching.aerial, coefs, tops, lgnd.dat, management_sf) {
  strt <- c(145.9,-14.5)
  survey_dates <- make_nice_date_range(bleacing.aerial) |>
    mutate(Region = factor(Region, levels = c("Northern", "Central", "Southern"))) |>
    arrange(Region)
  SurveyDates <- paste0("Survey dates\n(", survey_dates |> pull(survey_date), ")")

  levels(bleaching.aerial$CAT) <-  c('0%','>0% - 10%', '>10% - 30%', '>30% - 60%', '>60%')

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
make_oz_inset <- function(qld_sf, management_sf) {
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
  g <- g + ggspatial::annotation_scale(location = "bl",
    pad_x = unit(0.02, "npc"),
    pad_y = unit(0.05, "npc"),
    width_hint = 0.15)
  g <- g + ggspatial::annotation_north_arrow(location = "bl",
    height = unit(0.15, "npc"),
    width = unit(0.1, "npc"),
    pad_x = unit(0.02 + (0.15 / 2) - (0.1 / 2) - 0.005, "npc"),
    pad_y = unit(0.1, "npc"),
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering)
g
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
                                 y=unit(0.65, 'npc'),
                                 ## y=unit(im.y, 'npc'),
                                 vjust = 1,
                                 hjust = 0.5,
                                 width=unit(im.width, 'npc')))
                                 ## width=unit(0.1, 'npc')))
}



#######################################################################
## The following function defines the base plot (which includes      ##
## coral cover)                                                      ##
## Inputs:                                                           ##
##   - manta.mod: a dataframe containing the annual coral cover and  ##
##                change for each reef                               ##
## Returns                                                           ##
##   a ggplot object representing the map base                       ##
#######################################################################
make_base_plot <- function(manta.mod,
                           gbr_sf,
                           qld_sf,
                           management_sf,
                           tops,
                           lgnd.dat,
                           towns12,
                           Zones,
                           SurveyDates,
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
  g_base <-
    manta.mod |>
    filter(REPORT_YEAR == final_year) |>
    mutate(cCover = cut(Cover, breaks = c(0, 0.10, 0.30, 0.50, 0.75, 1), labels = 1:5)) |>
    ggplot(aes(y = Latitude, x = Longitude)) +
    geom_sf(data = gbr_sf |> filter(FEAT_NAME == "Reef"),
      fill = "grey",
      colour = "grey70", inherit.aes = FALSE) +
    geom_sf(data = qld_sf, fill = "grey", colour = "grey70", inherit.aes = FALSE) +
    geom_segment(data = grd.y, aes(x = x, y = y, xend = xend, yend = yend),
      size = 0.1, color = "grey80") +
    geom_segment(data = grd.x, aes(x = x, y = y, xend = xend, yend = yend),
      size = 0.1, color = "grey80") +
    annotate(geom = "text", x = -Inf, y = seq(-23, -11, by = 2),
      label = degreeLabelsNS(seq(-23, -11, by = 2)), hjust = -0.1, parse = TRUE, size = 3) +
    annotate(geom = "text", y = -Inf, x = seq(142, 154, by = 2),
      label = degreeLabelsEW(seq(142, 154, by = 2)), vjust = -0.3, parse = TRUE, size = 3) +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.5,
      inherit.aes = FALSE) +
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
    ## ## Northern, Central, Southern labels
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
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  g_base
}

## ---- Change plot
make_change_plot <- function(manta.mod, management_sf, tops, lgnd.dat, coefs) {
  strt <- c(145.9,-14.5)
  manta_mod_reefs <- manta.mod |>
    ungroup() |>
    filter(REPORT_YEAR == final_year) |>
    distinct() |>
    dplyr::select(REEF_NAME)
  lgnd.dat = data.frame(
    x = c(145.2, 145.2, 145.2),
    y = c(-11.3, -12, -12.7)
  )
  lgnd.dat$lab <- c("Increase","No~change","Decrease")
  lgnd.dat$Cat <- factor(c("increase","neutral", "decrease"))
  ## lgnd.dat.size <- ff(lgnd.dat[2,1:2],coefs,length<-5)
  lgnd.dat.size <- ff(c(coefs[1] + ((strt[2]+0.7)*coefs[2]), (strt[2]+0.7)),coefs,length<-5)
  lgnd.dat.size$lab <- format(seq(2.5,12.5,by<-2.5),digits<-2)
  lgnd.dat.size$size <- seq(2.5,12.5,by<-2.5)

  ly.mod <- manta.mod |>
    mutate(
      Diff = Median,
      D = D
    ) |>
    right_join(manta_mod_reefs) |>
    filter(!is.na(REPORT_YEAR))
  ## ly.mod |> write.csv(file='../data/processed/FinalYear_coral_change_2023.csv')

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
    filter(variable == "AnnualDiff") |>
    ## mutate(D = factor(D, levels = c("increase", "neutral", "decrease"))) |>
    mutate(D = factor(D, levels = c("decrease", "no change", "increase"))) |>
    ggplot(aes(y = Latitude, x = Longitude)) +
    geom_sf(data = management_sf, fill = NA, colour = "black", size = 0.2, inherit.aes = FALSE) +
    geom_point(aes(y = Latitude,x = Longitude,fill = D, size = abs(Diff*100)),
      alpha = 1,shape = 21,color = 'black', show.legend  =  TRUE)  +
    scale_fill_manual('Change', breaks = c("increase", "neutral", "decrease"),
      labels = c('Increase',"No change", 'Decrease'),
      values = c('green',"gray", 'red'),limits = c("increase", "neutral", "decrease")) +
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
  tem <- getwd()
  setwd(tempdir())
  g.change <- ggplotGrob(g.change)
  setwd(tem)
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

CoralTrends_generate_compilation_figure <- function(coral_cover_and_change, cots_and_bleaching) {
  require(sf)
  ## qld_sf <- st_as_sf(qld)
  gbr_sf <- st_read(dsn = "../data/spatial/Features/Great_Barrier_Reef_Features.shp")
  qld_sf <- get(load(file='../data/spatial/qld.RData')) |>
    st_as_sf()

  components_list <- list()

  coral_cover_and_change <- coral_cover_and_change |>
    CoralTrends_process_manta_mod_data__()

  spatial_3Zone <- readRDS(file = paste0(spatial_path, "spatial_3Zone.rds"))
  ## gbr_3Zone <- readRDS(file = paste0(spatial_path, "gbr_3Zone.rds"))
  management_sf <- spatial_3Zone |> st_as_sf()
  management_df <- get_gbr_df(spatial_3Zone |> st_as_sf())
  tops <- get_gbr_tops(management_df)
  coefs <- get_management_coefs(management_df)
  lgnd.dat <- get_lgnd_pos(coefs)

  if (all(is.na(cots_and_bleaching$Aerial_bleaching))) {  ## E.g. no aerial bleaching this year
    extra_width <- 18
  } else {
    extra_width <- 23
  }

  ## 1. Base plot
  survey_dates <- cots_and_bleaching |>
    mutate(Date = as.Date(`in-water_sample_date`, format = "%d/%m/%Y")) |>
    mutate(Region = factor(Region, levels = c("Northern GBR", "Central GBR", "Southern GBR"))) |>
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

  SurveyDates <- paste0("Survey dates\n(", survey_dates |> pull(survey_date), ")")
  manta.mod <- coral_cover_and_change

  g_base <- make_base_plot(
    manta.mod |> filter(variable == "Cover"),
    gbr_sf,
    qld_sf,
    management_sf,
    tops,
    lgnd.dat,
    towns12,
    Zones = c("Northern", "Central", "Southern"),
    SurveyDates,
    extra_width = extra_width
  )
  components_list$coral_cover <- manta.mod |>
    filter(variable == "Cover") |>
    dplyr::select(Latitude, Longitude, Cover) |>
    mutate(cCover = cut(Cover, breaks = c(0, 0.10, 0.30, 0.50, 0.75, 1), labels = 1:5))

  ## 2. Coral change
  g_change <- make_change_plot(
    manta.mod |> filter(variable == "AnnualDiff"),
    management_sf,
    tops,
    lgnd.dat,
    coefs)

  g <- add_change_to_plot(g = g_base,
    g_change,
    management_sf,
    shiftacross = 5.6)

  components_list$coral_cover_change <- manta.mod |>
    filter(variable == "AnnualDiff") |>
    dplyr::select(Latitude, Longitude, Diff = Median, D)

  ## 3. COTS
  cots <- process_cots(cots_and_bleaching)
  g_cots <- make_cots_plot(cots, coefs, tops, lgnd.dat, management_sf)
  g <- add_change_to_plot(g = g,
    g_cots,
    management_sf,
    shiftacross = 11.2)
  components_list$cots <- cots |>
    dplyr::select(
      Latitude = REEF_LAT, Longitude = REEF_LONG,
      OUTBREAK.CAT
    )

  ## 4. Bleaching
  bleaching <- process_bleaching(cots_and_bleaching, cots)
  g_bleaching <- make_bleaching_plot(bleaching, coefs, tops, lgnd.dat, management_sf)
  g <- add_change_to_plot(g = g, g_bleaching, management_sf, shiftacross = 16.8)
  components_list$bleaching <- bleaching |>
    dplyr::select(
      Latitude = REEF_LAT, Longitude = REEF_LONG,
      CAT = CAT
    )

  ## ## 4b. Aerial Bleaching
  ## bleaching_aerial <- process_aerial_bleaching(cots_and_bleaching)
  ## g_bleaching_aerial <- make_aerial_bleaching_plot(bleaching_aerial, coefs, tops, lgnd.dat, management_sf)
  ## g <- add_change_to_plot(g = g, g_bleaching_aerial, management_sf, shiftacross = 22.4)
  ## components_list$aerial_bleaching <- bleaching_aerial |>
  ##   dplyr::select(
  ##     Latitude = REEF_LAT, Longitude = REEF_LONG,
  ##     CAT
  ##   )

  ## 5. Add the Australia inset
  g_oz <- make_oz_inset(qld_sf, management_sf)
  g <- add_oz_to_plot(g, g_oz, management_sf, extra_width = extra_width,
    inset_width = 5, xmin = 170, ymin = -15)
  ## 6. Add the scalebar and north arrow
  g <- add_scale_and_arrow_to_plot(g)
  ## 7. AIMS watermark
  water_mark <- CoralTrends_get_aims_logo(wch = "AIMS stacked")
  g <- add_watermark_to_plot(g, water_mark, management_sf,
    extra_width = extra_width, oz_width = 5, inset_width = 3)

  ## 8. Save component_list
  ## saveRDS(components_list, file = '../data/modelled/map_components-2025.rds')


  fname <- paste0(fig_path, "cover_change_bleaching_cots.png")
  ggsave(filename = fname, plot = g, width = 15, height = 15 / 1.7)
  ggsave(filename = gsub(".png", ".hires.png", fname), plot = g, width = 15, height = 15 / 1.7, dpi = 300)
  ggsave(filename = gsub(".png", ".pdf", fname), plot = g, device = cairo_pdf, width = 15, height = 15 / 1.7)
  zip(zipfile = gsub(".png", ".zip", fname),
    files = c(fname, gsub(".png", ".hires.png", fname), gsub(".png", ".pdf", fname)),
    flags = "-j"
  )
}
