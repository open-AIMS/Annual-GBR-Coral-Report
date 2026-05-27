## To run interactively (for debugging)
## module load singularity
## singularity exec -B /export/project/monitoring-dashboard/Annual-GBR-Coral-Report:/home/Project -B /etc/localtime:/etc/localtime -B /etc/timezone:/etc/timezone -B /export/project/monitoring-dashboard/data:/data -B /export/project/monitoring-dashboard/output:/output --pwd /home/Project/R annual_report.sif R
## setwd("/export/project/monitoring-dashboard/Annual-GBR-Coral-Report/")
## Fake the args
if (1 == 2) { #for debugging only
  args <- c("--data_path=/data/",
    "--output_path=/output/",
    "--path=/data",
    "--path='/data/annual_report/2021-01-14/process/ALL/2024/Northern GBR/Regions/Northern GBR/raw/reef_data.csv'",
    "--method=annual_report",
    "--domain=Northern GBR",
    "--scale=region",
    "--refit=TRUE"
  )
  args <- c("--data_path=/data/",
    "--output_path=/output/",
    "--path='/data/annual_report/2021-01-14/process/ALL/2024/GBR/Regions/GBR/raw/reef_data.csv'",
    "--method=annual_report",
    "--domain=GBR",
    "--scale=region",
    "--refit=TRUE"
  )
}
setwd("../R")

args <- commandArgs()

cat(args)
cat("\n", getwd(), "\n")
cat("\nFiles in R:\n",
  paste("\t-", list.files(), collapse = "\n"),
  "\n")

cat("\nFiles in ../data\n",
  paste("\t-", list.files("../data"), collapse = "\n"),
  "\n")

cat("\nFiles in /data\n",
  paste("\t-", list.files("/data"), collapse = "\n"),
  "\n")

#########################################################################
## Read in the project specific functions                              ##
## This also includes running a function (CoralTrends_checkPackages()) ##
## that assesses whether all the required packages are present         ##
#########################################################################
source('CoralTrends_functions.R')
CoralTrends_checkPackages()

CoralTrends_parse_args(args)

#########################################################################
## Read and load the configurations set in parameters/CoralTrends.conf ##
#########################################################################
source('CoralTrends_config.R')

######################################################################
## Generate the three zone configuration based on De'ath et al 2012 ##
##    Zones are stored in:                                          ##
##    - data/spatial/whagbr.RData                                   ##
##    - data/spatial/whagbr.n.RData                                 ##
##    - data/spatial/whagbr.c.RData                                 ##
##    - data/spatial/whagbr.s.RData                                 ##
##    - data/spatial/qld.RData                                      ##
######################################################################
source('CoralTrends_spatial_3Zone.R')

#######################################################################
## Retrieve the data from the path provided by the command line      ##
## arguments.  This replaces the old function that extracted the data##
## from the oracle database.
#######################################################################
## source('CoralTrends_getData_Manta.R')

if (file.exists(path)) {
  data <- read_csv(gsub(".zip", ".csv", path), trim_ws = TRUE)
  saveRDS(data, file = paste(data_path, "data_", domain, ".rds"))
} else stop("Raw data does not exist")

#################################################
## Process the Manta-tow data                  ##
## - Only include data collected after 1985    ##
## - Convert data into percent cover           ##
## - Summarize the data to the reef/year level ##
## - Assign a zone                             ##
## Data are stored in:                         ##
## - data/processed/manta.sum.RData            ##
#################################################
source('CoralTrends_processData_Manta_3Zone.R')


###############################
## Fit STAN (or INLA) models ##
###############################
source('CoralTrends_stan_Manta_3Zone.R')

######################################################################
## Generate single panel coral cover trend plot                     ##
## ggplot object is stored in:                                      ##
## - data/modelled/gg_data_domain_COVER_beta_ _ _ _Cover.rds        ##
## ggplot figure is stored in:                                      ##
## - data/modelled/gg_fig_domain_COVER_beta_ _ _ _Cover.png         ##
######################################################################
source('CoralTrends_trend_Manta_3Zone.R')





cat("\n\n\n", append = TRUE)
