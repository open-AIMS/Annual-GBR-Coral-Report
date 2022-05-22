source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

load("../data/spatial/whagbr.RData")
######################################################################
## Clip the world heritage area gbr object to latitudinal divisions ##
## defined by De'ath 2012                                           ##
######################################################################
whagbr.n = ML_gClip(whagbr, cbind(c(142,-15.4),c(155,0)))
#whagbr.c = ML_gClip(whagbr, cbind(c(142,-20),c(155,-15.4)))
bb=rbind(c(142,-20.7),
         c(148.7,-20.7),
         c(152,-19.6),
         c(152,-15.4),
         c(142,-15.4))
b.poly=SpatialPolygons(list(Polygons(list(Polygon(bb)),ID=1)))
whagbr.c = ML_gClip(whagbr, b.poly)
#whagbr.s = ML_gClip(whagbr, cbind(c(142,-25),c(155,-20)))
bb=rbind(c(142,-20.7),
         c(148.7,-20.7),
         c(152,-19.6),
         c(155,-19.6),
         c(155,-25),
         c(142,-25))
b.poly=SpatialPolygons(list(Polygons(list(Polygon(bb)),ID=1)))
whagbr.s = ML_gClip(whagbr, b.poly)

## data(qld)
load("../data/spatial/qld.RData")

## save(whagbr, file='../data/spatial/whagbr.RData')
save(whagbr.n, file='../data/spatial/whagbr.n.RData')
save(whagbr.c, file='../data/spatial/whagbr.c.RData')
save(whagbr.s, file='../data/spatial/whagbr.s.RData')
## save(qld, file='../data/spatial/qld.RData')

##Consolidate into one spatial_3Zone
whagbr.n=spChFIDs(whagbr.n,'N')
whagbr.c=spChFIDs(whagbr.c,'C')
whagbr.s=spChFIDs(whagbr.s,'S')
spatial_3Zone = rbind(whagbr.n, whagbr.c, whagbr.s)
#spatial_3Zone = SpatialPointsDataFrame(spatial_3Zone
save(spatial_3Zone, file='../data/spatial/spatial_3Zone.RData')
