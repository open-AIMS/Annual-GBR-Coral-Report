########################################################################
## This module loads a cofiguration file that contains paramaters for ##
## settings to be applied across all aspects of the analyses          ##
## This file is a text file with key=value pairs                      ##
##                                                                    ##
## Specifically, the pairs are:                                       ##
## Size:                    the number of bootstrapp samples          ##
## Timeseries.plot.width    the width of a timeseries plot            ##
## EndDate                  the last day (date) of the current focal  ##
##                            year (must be in YYYY-MM-DD format)     ##
## Timeseries.plot.height   the height of a timeseries plot           ##
## FocalYear                the current report card year              ##
## StartDate                the lower date range cutoff (must be in   ##
##                            YYY-MM-DD format)                       ##
########################################################################


CoralTrends_tryCatch(
    {
        if(!dir.exists('../data')) dir.create('../data')
        if(!dir.exists('../data/primary')) dir.create('../data/primary')
        if(!dir.exists('../data/processed')) dir.create('../data/processed')
        if(!dir.exists('../data/spatial')) dir.create('../data/spatial')
        if(!dir.exists('../data/modelled')) dir.create('../data/modelled')
        
        if(!dir.exists('../output')) dir.create('../output')
        if(!dir.exists('../output/figures')) dir.create('../output/figures')
        if(!dir.exists('../output/tables')) dir.create('../output/tables')
        
        if(!dir.exists('../logs')) dir.create('../logs')
 
        return=NULL
    }, '../logs/all.log','--Config--',msg='configure necessary folders', return=NULL)

## ---- Config
CoralTrends_tryCatch(
    {
        config = readLines('../parameters/CoralTrends.conf')
        config = gsub('(.*Date)=(.*)','\\1=as.Date(\'\\2\')',config)
        eval(parse(text=config))
        return=NULL
    }, '../logs/all.log','--Config--',msg='load general configurations', return=NULL)
## ----

