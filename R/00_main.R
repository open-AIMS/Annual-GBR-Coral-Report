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






cat("\n\n\n", append = TRUE)
