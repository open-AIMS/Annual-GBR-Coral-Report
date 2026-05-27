cat(paste0("\nGenerate nice trend for ",
  domain, " data\n====================================\n"),
  append = TRUE)

source('CoralTrends_functions.R')
source('CoralTrends_trend_functions.R')

library(sf)
library(oz)

data <- readRDS(file = paste0(data_path, "modelled/annual_report_region_", domain,
  ".rds"))

## gbr_3Zone <- readRDS(file = paste0(spatial_path, "gbr_3Zone.rds"))

hues <- RColorBrewer::brewer.pal(4, "Blues")

## Generate mimimap banner

## cat(paste0("\t - generate the minimap banner\n"),
##   append = TRUE)
## data <- data |> CoralTrends_generate_single_banner(region = domain)
## ggsave("test.png",
##   data$minimap[[1]] + theme(plot.background = element_rect(fill = "white")),
##   width = 6, height = 6)

## ## Make the trend plot
cat(paste0("\t - generate the single trend plot\n"),
  append = TRUE)
data <- data |> CoralTrends_generate_single_trend(
  region = domain,
  final_year = NA
)
## ## ggsave("test.png",
## ##   plot = trend_plot,
## ##   width = 6, height = 6)

## ## Finalise plot (add banner to trend plot)
cat(paste0("\t - finalise the plot (add minimap to facet)\n"),
  append = TRUE)
data <- data |>
  CoralTrends_addorn_plots(region = domain)




fname <- paste0(data_path, "modelled/annual_report_", domain, "_with_figs.rds")
saveRDS(data, fname)

cat(paste0("\t - Finished!\n"),
  append = TRUE)
