cat("\nInside CoralTrends_compile_figures.R\n=============================================\n")

cat("The current working directory is: ", getwd(), "\n\n")

args <- commandArgs()

## Fake args for testing
if (length(args) == 1)  { ## clearly being run interactively
  args <- c(
    "--cwd=../../Annual-GBR-Coral-Report/R/",
    "--domain=annual_report",
    "--data_path=/export/project/monitoring-dashboard/data/",
    "--final_year=2026",
    "--output_path=/output/",
    "--sql_path=../data/",
    "--path=/data/",
    "--refit=TRUE",
    "--log=/export/project/monitoring-dashboard/data/log/"
  )
}

cat(args)

has_cwd_argument <- any(grepl("--cwd=.*", args, perl = TRUE))
if(has_cwd_argument) {
  arg <- args[grep("--cwd=.*", args)]
  cwd <- gsub("--cwd=(.)", "\\1", arg)
  setwd(cwd)
} else stop("no --cwd= component provided")
cat("\nNew working directory is: ", getwd(), "\n")

## #########################################################################
## ## Read in the project specific functions                              ##
## ## This also includes running a function (CoralTrends_checkPackages()) ##
## ## that assesses whether all the required packages are present         ##
## #########################################################################
source('CoralTrends_functions.R')
source('CoralTrends_compile_figures_functions.R')
source('../../dev/R/backend_functions.R')

CoralTrends_checkPackages()
library(magick)
library(png)
library(DBI)

CoralTrends_parse_args(args)


cat("\n\n", ls())

config_ <- list(
  dashboard_log = log_file,
  db_path = "/export/project/monitoring-dashboard/data/dashboard.sqlite",
  model_path = "/export/project/monitoring-dashboard/data/modelled/"
  )


########################################################
## get the raw data compiliation
cat("\n- get the raw data ", "\n")
raw_manta <- CoralTrends_get_raw_data_via_db_and_file() |>
  CoralTrends_tidy_raw_data()

cat("\n- get the annual estimate summaries data ", "\n")
annual_manta <- CoralTrends_get_annual_data_via_db_and_file(type = "year_sum") |>
  CoralTrends_tidy_annual_manta_data(raw_data, type = "year_sum")

cat("\n- get the annual contrast estimate summaries data ", "\n")
annual_yearcomp_manta <- CoralTrends_get_annual_data_via_db_and_file(type = "yearcomp_sum") |>
  CoralTrends_tidy_annual_manta_data(raw_data, final_year = final_year, type = "yearcomp_sum")

cat("\n- combine cover and change", "\n")
coral_cover_and_change <- CoralTrends_combine_cover_and_change(annual_manta,
  annual_yearcomp_manta, final_year = final_year)

cat("\n- coral_cover_and_change",capture.output(print(head(coral_cover_and_change))), "\n", sep = "\n")

## sql_path <- "../data/"
cat("\n- extract COTS and bleaching", "\n")
cots_and_bleaching <- CoralTrends_get_cots_and_bleaching()


cat("\n- generate compilation figure", "\n")
CoralTrends_generate_compilation_figure(coral_cover_and_change, cots_and_bleaching)


cat("\nEnd of CoralTrends_compile_figures.R\n=============================================\n")
