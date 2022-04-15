require(sf)
library(tidyverse)
library(raster)
library(tmap)

### extent of region ####

#EPSG=32632
Denmark <- getData('GADM', country='DNK',level=0) %>%
            st_as_sf() %>%
            st_transform(.,crs=32632)

### full grid ####

newres = 1000#1 km grid
mygrid <- raster(extent(as(Denmark,'Spatial')))
res(mygrid) <- newres 
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
plot(Denmark,add=T)

### mask grid #####

mygrid <- mask(mygrid, Denmark)
plot(mygrid)

### take samples ####

mypolys <- raster::rasterToPolygons(mygrid)
mypoints<- raster::rasterToPoints(mygrid)

#get land use data
ldir <- "/Users/dianabowler/Dropbox/DOF/landuseData/basemap03_2018"

#aggregated layer with all classes
r <- raster(paste(ldir,"lu_agg_2018.tif",sep="/"))
res(r)
r <- aggregate(r, fact=10,fun="modal")
r[r==999999] <- NA
#sort(unique(values(r)))#corresponds with table 3.6!!
#table(values(r))
plot(r)
crs(r) <- "epsg:32632"
plot(Denmark,add=T)

#overlay land map with squares
tm_shape(r)+
  tm_raster()+
  tm_shape(Denmark)+
  tm_borders()
#overlap fine

#extract land use for the kvatradnr
squares_buffer_landUse <- raster::extract(r,mypolys,df=T,weights=TRUE)
names(squares_buffer_landUse)[2] <- "LandUse"
squares_buffer_landUse <- subset(squares_buffer_landUse, !is.na(LandUse))

#summarize per ID
squares_buffer_landUse_summary <- squares_buffer_landUse %>%
  dplyr::group_by(ID,LandUse) %>%
  dplyr::summarise(weight=sum(weight))


#see Table 3.6 in TR159.pdf for explanation of land use codes
table(squares_buffer_landUse_summary$LandUse)

#now get proportion of urban areas and agricultural areas
landuse_summary <- squares_buffer_landUse_summary %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(urban = sum(weight[LandUse %in% 110000:125110]),
                   open_urban = sum(weight[LandUse %in% c(126000:130110,150000:160000)]),
                   agri_int = sum(weight[LandUse %in% 211000:212000]),
                   agri_ext = sum(weight[LandUse %in% c(220000,230000,321000,321220)]),
                   forest = sum(weight[LandUse %in% 311000:312000]),
                   wet = sum(weight[LandUse %in% 322000:322220]),
                   freshwater = sum(weight[LandUse %in% 411000:412000]),
                   seawater = sum(weight[LandUse %in% 420000]),
                   mapped = sum(weight[!LandUse %in% 800000]),
                   total = sum(weight))

head(landuse_summary)

#combine with land use data
squares_buffer <- bind_cols(mypolys@data,landuse_summary)
head(squares_buffer)
names(squares_buffer)[1] <- "grid"

#add on coordinates
squares_buffer$X <- coordinates(mypolys)[,1]
squares_buffer$Y <- coordinates(mypolys)[,2]

#plot to check
ggplot(squares_buffer)+
  geom_point(aes(x=X,y=Y,colour=urban))          

#also distance to coast
Germany <- getData('GADM', country='DEU',level=0) %>%
  st_as_sf() %>%
  st_transform(.,crs=32632)

coast <- st_union(bind_rows(Denmark,Germany))
plot(coast)

#get points
mygridPoints <- sf::st_as_sf(as.data.frame(mypoints), coords = c("x", "y"), 
                             crs = st_crs(Denmark)) 
ggplot()+
  geom_sf(data=coast)+
  geom_sf(data=mygridPoints,color="red")

#get distance between them
mygridPoints1 <- sf::st_transform(mygridPoints, crs=4326) %>%
                    sf::st_coordinates(mygridPoints)
coast1 <- sf::st_transform(coast, crs=4326) 
dist <- geosphere::dist2Line(p = as.matrix(mygridPoints1), 
                             line = as(coast1,'Spatial'))

#pack into data frame
dist <- as.data.frame(dist)
dist$grid <- mygridPoints$layer
names(dist)[1] <- "distCoast"

saveRDS(squares_buffer,file="environ-data/grid_environdata.rds")
