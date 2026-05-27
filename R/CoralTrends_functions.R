##################################################################
## The following code forms part of the Annual GBR Coral trends ##
## analysis code base.                                          ##
##################################################################


###################################################################
## The following function checks to ensure that all the required ##
## packages are available on the system.                         ##
###################################################################
CoralTrends_checkPackages <- function() {
    require(gdata) # load this first as it masks many
    require(tidyverse)
    require(gtable)
    require(grid)
    require(gridExtra)
    require(xtable)

    require(rgdal)
    require(rgeos)
   ## require(mapping) # consider replacing this with a self contained function
    require(maptools)
    require(raster)

    require(rstanarm)
    require(coda)
    require(glmmTMB)
    require(INLA)
    require(emmeans)
    library(mgcv)
    require(betareg)
}


#########################################################################
## The following function parses the command line arguments.           ##
## Note.  This is just about the only function that does not trigger a ##
## tryCatch log file - and this is primarily because before running    ##
## this function, the location of writable directories is unknown      ##
## since they are parsed as command line arguments.                    ##
#########################################################################
CoralTrends_parse_args <- function(args) {
  cat("\nParsed command line arguments\n==================================\n")
  ## external directories (parsed from command line)
  has_data_path_argument <- any(grepl("--data_path=.*", args, perl = TRUE))
  if(has_data_path_argument) {
    arg <- args[grep("--data_path=.*", args)]
    data_path <- gsub("--data_path=(.)", "\\1", arg)
    cat(paste0("\t- data path: ", data_path), "\n",
      append = TRUE)
    assign("data_path", data_path, envir = .GlobalEnv)
  } else stop("no --data_path= component provided")

  spatial_path <- paste0(data_path, "spatial/")
  if(!dir.exists(spatial_path)) dir.create(spatial_path)
  assign("spatial_path", spatial_path, envir = .GlobalEnv)

  has_domain_argument <- any(grepl("--domain=.*", args, perl = TRUE))
  if(has_domain_argument) {
    arg <- args[grep("--domain=.*", args)]
    domain <- gsub("--domain=(.)", "\\1", arg)
    cat(paste0("\t- domain: ", domain), "\n",
      append = TRUE)
    assign("domain", domain, envir = .GlobalEnv)
  } else stop("no --domain= component provided")

  log_path <- paste0(data_path, "log/")
  if(!dir.exists(log_path)) dir.create(log_path)
  log_file <- paste0(log_path, "annual_report_", domain, ".log")
  if (file.exists(log_file)) file.remove(log_file)
  assign("log_file", log_file, envir = .GlobalEnv)
    cat(paste0("\t- temporary log file: ", log_file), "\n",
      append = TRUE)

  has_output_path_argument <- any(grepl("--output_path=.*", args, perl = TRUE))
  if(has_output_path_argument) {
    arg <- args[grep("--output_path=.*", args)]
    output_path <- gsub("--output_path=(.)", "\\1", arg)
    cat(paste0("\t- output path: ", output_path), "\n",
      append = TRUE)
    assign("output_path", output_path, envir = .GlobalEnv)

    fig_path <- paste0(output_path, "figures/")
    cat(paste0("\t- figures path: ", fig_path), "\n",
      append = TRUE)
    assign("fig_path", fig_path, envir = .GlobalEnv)

    tbl_path <- paste0(output_path, "tables/")
    cat(paste0("\t- tables path: ", tbl_path), "\n",
      append = TRUE)
    assign("tbl_path", tbl_path, envir = .GlobalEnv)
  } else stop("no --output_path= component provided")

  has_path_argument <- any(grepl("--path=.*", args, perl = TRUE))
  if(has_path_argument) {
    arg <- args[grep("--path=.*", args)]
    path <- gsub("--path=(.)", "\\1", arg)
    path <- gsub("'", "", path)
    cat(paste0("\t- path: ", path), "\n",
      append = TRUE)
    assign("path", path, envir = .GlobalEnv)
  } else stop("no --path= component provided")

  has_refit_argument <- any(grepl("--refit=.*", args, perl = TRUE))
  if(has_refit_argument) {
    arg <- args[grep("--refit=.*", args)]
    refit <- gsub("--refit=(.)", "\\1", arg)
    refit <- ifelse(refit == "TRUE", TRUE, FALSE)
    cat(paste0("\t- refit: ", refit, " class:", class(refit)), "\n",
      append = TRUE)
    assign("refit", refit, envir = .GlobalEnv)
  } else stop("no --refit= component provided")

  has_final_year_argument <- any(grepl("--final_year=.*", args, perl = TRUE))
  if(has_final_year_argument) {
    arg <- args[grep("--final_year=.*", args)]
    final_year <- as.numeric(as.character(gsub("--final_year=(.)", "\\1", arg)))
    cat(paste0("\t- final_year: ", final_year), "\n",
      append = TRUE)
    assign("final_year", final_year, envir = .GlobalEnv)
  } else stop("no --final_year= component provided")

  has_sql_path_argument <- any(grepl("--sql_path=.*", args, perl = TRUE))
  if(has_sql_path_argument) {
    arg <- args[grep("--sql_path=.*", args)]
    sql_path <- gsub("--sql_path=(.)", "\\1", arg)
    cat(paste0("\t- sql_path: ", sql_path), "\n",
      append = TRUE)
    assign("sql_path", sql_path, envir = .GlobalEnv)
  }
}

#########################################################################
## The following function appends a log statement into a log file      ##
## parameters:                                                         ##
##     status:    a string indicating either 'FAILURE', 'SUCCESS' or   ##
##                'WARNING'                                            ##
##     logFile:   a character string representation of the log file    ##
##                name (including path relative to current working     ##
##                directory)                                           ##
##     Category:  a character string with a category to appear         ##
##                verbatim in the log                                  ##
##     msg1:      a character string with a message to appear verbatim ##
##                in the log                                           ##
#########################################################################
CoralTrends_log <- function (status, logFile='data/logs/env.log',Category, msg1) {
    options(digits.secs=2)              ## switch to subsecond display
    ## Check if the log file exists, and if it does not, create it
    d=dirname(logFile)
    files <- list.files(d)
    if(!any(grepl(paste0('^',logFile,'$'),files))) system(paste0('touch ',logFile))
    now <- Sys.time()

    msg <- paste0(now, '|',status, ': ', Category, ' ',msg1)
    if( !is.null(msg)){ write(msg,file=paste0(logFile),append=TRUE)}

}

##########################################################################
## The following function provides a more useful error handling         ##
## routine.                                                             ##
##    expr:      an R expression to be evaluated                        ##
##    logFile:   a character string represetnation of the log file name ##
##               (including path relative to the current working        ##
##               directory)                                             ##
##    Category:  a character string representation of error category    ##
##    msg:       a character string with a message to appear verbatim   ##
##               in the log                                             ##
##    return:    boolean, whether to return a TRUE or FALSE             ##
##########################################################################
CoralTrends_tryCatch <- function(expr, logFile,Category, expectedClass=NULL, msg=NULL, return=NULL, showWarnings=FALSE) {
    #msg <- paste0(now, '| ', msg)
    max.warnings<-4
    warnings<-0
    W <- NULL
    w.handler <- function(w){ # warning handler
        m<-w$message
        if ((warnings < max.warnings) && (grepl ('CoralTrends_WARNING', m)>0)) {
            CoralTrends_log('WARNING', logFile,Category, paste(warnings, msg, m))
            warnings<<-warnings+1
        }
        invokeRestart("muffleWarning")
    }
    ret <- list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                    warning = w.handler),warning = W)
    if(!is.atomic(ret$value) && !is.null(ret$value$message)){
        ## An error occurred
        class(ret) <- "try-error"
        CoralTrends_log('WARNING', logFile,Category, paste(msg, ret$value$message))
        if(!is.null(return)) {
            FALSE
        }#else return()
    } else {    #no error check for warning
        CoralTrends_log('INFO', logFile, Category, msg)
        if(!is.null(return)) {
            TRUE
        }
    }
}

#####################################################################
## The following function clips the a SpatialPolygons object to    ##
## particular latitudinal divisions so as to make the three De'ath ##
## 2012 zones.                                                     ##
##  parameters:                                                    ##
##    shp:  a SpatialPolygons object                               ##
##    bb:   a matrix of bounding box longitude and latitudes       ##
##  returns: a SpatialPolygons object                              ##
#####################################################################
ML_gClip <- function(shp, bb){
  if(any(class(bb) == "matrix")) {
      if (identical(dim(bb), c(2L,2L))) {
          b_poly <- as(raster:::extent(as.vector(t(bb))), "SpatialPolygons")
      } else b_poly = bb
  } else if (class(bb) =='SpatialPolygons') {
      b_poly=bb
  } else b_poly <- as(raster:::extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

#######################################################################
## The following function converts the Manta tow coral cover classes ##
## into percent cover values in the range [0,1].                     ##
##   parameters:                                                     ##
##      x:     a character vector of coral class categories          ##
##   returns:  a numeric cover abundance                             ##
#######################################################################
CoralTrends_calcPercent = function(x) {
    ifelse(x=='0', 0,
    ifelse(x=='1', 0.05,
    ifelse(x=='1L', 0.025,
    ifelse(x=='1U', 0.075,
    ifelse(x=='2', 0.2,
    ifelse(x=='2L', 0.15,
    ifelse(x=='2U', 0.25,
    ifelse(x=='3', 0.4,
    ifelse(x=='3L', 0.35,
    ifelse(x=='3U', 0.45,
    ifelse(x=='4', 0.625,
    ifelse(x=='4L', 0.5625,
    ifelse(x=='4U', 0.6875,
    ifelse(x=='5', 0.875,
    ifelse(x=='5L',0.8125,0.9375)))))))))))))))
}

## The following function converts a cover into a category representing the median of the category
median_cover_cat <- function(dat) {
    n <- length(dat)
    if (n %% 2 == 0) {
        d <- paste(unique(sort(dat)[c((n-1)/2, (n+1)/2)]), collapse='/')
    } else {
        d <- paste(unique(sort(dat)[(n+1)/2]))
    }
    factor(d)
}



######################################################################
## The following function generates Zone categories based on De'ath ##
## 2012's latitudinal divisions.                                    ##
##   parameters:                                                    ##
##      x:     a numeric vector of latitudes                        ##
##   returns:  a categorical vector of Zones                        ##
######################################################################
CoralTrends_calc3ZoneLocations <- function(x) {
    factor(ifelse(x<= -10.68 & x > -15.4, 'Northern',  #glenn's version is -11.8
           ifelse(x<= -15.4 & x > -20.0, 'Central',
           ifelse(x<= -20.0 & x > -23.92, 'Southern','Outside'))))
}

CoralTrends_calc3ZoneLocation <- function(dat) {
  gbr_3Zone <- readRDS(file = paste0(spatial_path, "gbr_3Zone.rds"))
  dat |>
    sf::st_as_sf(coords = c("Longitude", "Latitude"),
      crs = sf::st_crs(gbr_3Zone)) |>
    sf::st_join(gbr_3Zone) %>%
    cbind(Longitude = sf::st_coordinates(.)[,1],
      Latitude = sf::st_coordinates(.)[,2]) |>
    sf::st_drop_geometry()
}

CoralTrends_get_aims_logo <- function(wch = "AIMS") {
  if (wch == "AIMS") {
    img <- magick::image_read(
      path = "../parameters/ECM_1280945_v1_AIMS Logo Stacked white 1200px (2).png"
    )
  } else if (wch == "AIMS stacked") {
    img <- magick::image_read(
      path = "../parameters/ECM_1280725_v1_AIMS Logo - stacked.jpg"
    )
  } else {
    img <- magick::image_read(
      path = "../parameters/AIMSLogo_Colour_inline.png"
    )
  }
  return(img)
}
