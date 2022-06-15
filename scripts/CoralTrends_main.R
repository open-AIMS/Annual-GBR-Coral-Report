## Please cross reference this against the codes written for the "Products for ALT" project

## The following code base generates the Annual GBR Coral trends 
## analysis.  This file constitutes the main script from which all
## others are called.  The code base was written by Murray Logan
## (AIMS) April 2018.
##
## BEFORE RUNNING THIS SCRIPT, PLEASE REVIEW:
## parameters/CoralTrends.conf
##
## The code base comprises the following:
## |
## |- CoralTrends_main.R
## |- CoralTrends_functions.R
## |- CoralTrends_config.R
## |- CoralTrends_spatial_3Zone.R
## |- CoralTrends_getData_Manta.R
## |- CoralTrends_processData_Manta_3Zone.R
## |- CoralTrends_spatialMap_3Zone.R
## |- CoralTrends_temporalHeatMap_3Zones.R
## |- CoralTrends_stan_Manta.R
## |- CoralTrends_trend_Manta_3Zone.R
## |- CoralTrends_map_Manta.R


#########################################################################
## Read in the project specific functions                              ##
## This also includes running a function (CoralTrends_checkPackages()) ##
## that assesses whether all the required packages are present         ##
#########################################################################
source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

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

######################################################################################################################################################################
## Extract and read in the Manta-tow data                                                                                                                           ##
## and store it under data/primary/manta.csv and data/primary/manta.RData                                                                                           ##
## The SQL is as follows                                                                                                                                            ##
## SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME,V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.REEF_LAT,      ##
##   V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL, V_RM_SAMPLE.SAMPLE_CLASS                                                                    ##
## FROM RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID                                                                               ##
## WHERE (((V_RM_SAMPLE.SAMPLE_CLASS) In ('K','C','G','Z') Or (V_RM_SAMPLE.SAMPLE_CLASS) Is Null))                                                                  ##
## ORDER BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR, ##
## RM_MANTA.TOW_SEQ_NO                                                                                                                                              ##
######################################################################################################################################################################
source('CoralTrends_getData_Manta.R')

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


###################################
## Generate site maps            ##
## Maps are stored in:           ##
## - output/figures/MapOfSites.* ##
###################################
source('CoralTrends_spatialMap_3Zone.R')

######################################################
## Generate temporal head maps of data availability ##
## Heat maps are stored in:                         ##
## - output/figures/TemporalHeatMap_3Zone.*        ##
######################################################
source('CoralTrends_temporalHeatMap_3Zone.R')

###############################
## Fit STAN (or INLA) models ##
###############################
## This needs to be run manually at the moment!
if (rerun_models) source('CoralTrends_stan_Manta_3Zone.R')

## 2022 - the rest is run locally
## will need a copy of the entire /data/modelled folder

############################
## Construct trend graphs ##
############################
## The following is the new version with just the wanted graphs
source('CoralTrends_summary_figures.R')
## The following makes all graphs and is getting dated
source('CoralTrends_trend_manta_3Zone.R')


#################################
## Construct coral change maps ##
#################################
source('CoralTrends_coralchange_map_3Zone.R')

source('CoralTrends_coralchange_more.R')


source('CoralTrends_recent_declines.R')


zip('output/data4Mike.zip',
    files=c(
        'output/figures/FourPlots1.pdf',
        'output/figures/FourPlots1.png',
        #'output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.pdf',
        #'output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.png',
        'output/figures/mantaCAT_northern.pdf',
        'output/figures/mantaCAT_northern.png',
        'output/figures/mantaCAT_central.pdf',
        'output/figures/mantaCAT_central.png',
        'output/figures/mantaCAT_southern.pdf',
        'output/figures/mantaCAT_southern.png',
        'output/tables/fullannualTable.RData',
        'data/modelled/ly.RData',
        'output/figures/3ZonesStanRibbons.pdf',
        'output/figures/3ZonesStanRibbons.png',
        'output/figures/3ZonesStanErrorBars.pdf',
        'output/figures/3ZonesStanErrorBars.png',
        'output/figures/3ZonesStan.pdf',
        'output/figures/3ZonesStan.png'),
    flags=c('-j','-o')
    )
