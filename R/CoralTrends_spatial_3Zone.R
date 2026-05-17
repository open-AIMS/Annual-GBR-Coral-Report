cat(paste0("\nDefining spatial domains ",
  domain, "data\n====================================\n"),
  append = TRUE)

source('CoralTrends_functions.R')
CoralTrends_checkPackages()

cat(paste0("\t- defining each region (zone)\n"),
  append = TRUE)
load("../data/spatial/whagbr.RData")
######################################################################
## Clip the world heritage area gbr object to latitudinal divisions ##
## defined by De'ath 2012                                           ##
######################################################################
whagbr.n <- ML_gClip(whagbr, cbind(c(142,-15.4),c(155,0)))
bb <- rbind(c(142,-20.7),
         c(148.7,-20.7),
         c(152,-19.6),
         c(152,-15.4),
         c(142,-15.4))
b.poly <- SpatialPolygons(list(Polygons(list(Polygon(bb)), ID = 1)))
whagbr.c <- ML_gClip(whagbr, b.poly)
bb <- rbind(c(142,-20.7),
         c(148.7,-20.7),
         c(152,-19.6),
         c(155,-19.6),
         c(155,-25),
         c(142,-25))
b.poly <- SpatialPolygons(list(Polygons(list(Polygon(bb)), ID = 1)))
whagbr.s <- ML_gClip(whagbr, b.poly)

load("../data/spatial/qld.RData")

cat(paste0("\t- exporting\n"),
  append = TRUE)
saveRDS(whagbr, file = paste0(spatial_path, "whagbr.rds"))
saveRDS(whagbr.n, file = paste0(spatial_path, "whagbr.n.rds"))
saveRDS(whagbr.c, file = paste0(spatial_path, "whagbr.c.rds"))
saveRDS(whagbr.s, file = paste0(spatial_path, "whagbr.s.rds"))

##Consolidate into one spatial_3Zone
whagbr.n <- spChFIDs(whagbr.n, "N")
whagbr.c <- spChFIDs(whagbr.c, "C")
whagbr.s <- spChFIDs(whagbr.s, "S")
spatial_3Zone <- rbind(whagbr.n, whagbr.c, whagbr.s)
saveRDS(spatial_3Zone, file = paste0(spatial_path, "spatial_3Zone.rds"))

cat(paste0("\t- creating gbr_3Zone sf object\n"),
  append = TRUE)
gbr_3Zone <- sf::st_as_sf(whagbr.n + whagbr.c + whagbr.s) %>%
  sf::st_set_crs('EPSG:4283') %>%
  mutate(Region=c('Northern GBR','Central GBR','Southern GBR'))
saveRDS(gbr_3Zone, file = paste0(spatial_path, "gbr_3Zone.rds"))
print(head(gbr_3Zone))
