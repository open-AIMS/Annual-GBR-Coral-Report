library(tidyverse)

## SQL for extracting manta tow data from oracle
## ---- Manta.SQL

## In 2023 there was a database issue that means that it was not
## considered safe to extract all the data from the database.
## Instead, it was decided that I should retain last years data and
## add only this years extraction on top of the past data.  When this
## database issue is addressed, the following conditional section
## should be removed.

writeLines("
SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME,V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG,
  V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL, V_RM_SAMPLE.SAMPLE_CLASS
FROM RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID
WHERE (((V_RM_SAMPLE.SAMPLE_CLASS) In ('K','C','G','Z') Or (V_RM_SAMPLE.SAMPLE_CLASS) Is Null))
ORDER BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LAT,
  V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO","../data/primary/manta.sql")
## ----

if (goto_database_manta) system("java -jar dbExport.jar ../data/primary/manta.sql ../data/primary/manta.csv reef reefmon") 

manta <- read.csv('../data/primary/manta.csv',strip.white=TRUE)  


## if (final_year == 2023) {

##     writeLines("
## SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME,V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG,
##   V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL, V_RM_SAMPLE.SAMPLE_CLASS
## FROM RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID
## WHERE (((V_RM_SAMPLE.SAMPLE_CLASS) In ('K','C','G','Z') Or (V_RM_SAMPLE.SAMPLE_CLASS) Is Null)) AND V_RM_SAMPLE.VISIT_NO LIKE '%31%'
## ORDER BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LAT,
##   V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO","../data/primary/manta_.sql")

##     if (goto_database_manta) system("java -jar dbExport.jar ../data/primary/manta_.sql ../data/primary/manta_.csv reef reefmon") 

##     manta_ <- read.csv('../data/primary/manta_.csv',strip.white=TRUE)  
##     manta <- manta %>% bind_rows(manta_)
## }

save(manta, file='../data/primary/manta.RData')

