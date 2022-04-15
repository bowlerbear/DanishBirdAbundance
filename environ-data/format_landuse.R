library(sf)
#library(raster)
library(terra)
library(tmap)
library(tidyverse)

### SQUARES ####

#read in dof squares
squares <- st_read(dsn = "data/Insects_and_TTT",
                   layer = "transect squares utm32")
plot(squares)

#group into kvadratnr
squares <- squares %>%
  dplyr::group_by(kvadratnr) %>%
  dplyr::summarise() %>%
  dplyr::ungroup() 

plot(squares)

#buffer the squares
summary(st_area(squares))#area is 1000000

#no buffer
squares_buffer <- squares

#5km buffer
#we want area is increase by 5, so length to increase by 2.236
squares_buffer <- st_buffer(squares,dist=659)

#check
summary(st_area(squares_buffer))


#get land use data
ldir <- "C:/Users/db40fysa/Dropbox/DOF/landuseData/basemap03_2018"
#aggregated layer with all classes
r <- rast(paste(ldir,"lu_agg_2018.tif",sep="/"))
res(r)
r <- aggregate(r, fact=10,fun="modal")
r[r==999999] <- NA
#sort(unique(values(r)))#corresponds with table 3.6!!
#table(values(r))
plot(r)
crs(r) <- "epsg:32632"


#overlay land map with squares
tm_shape(raster::raster(r))+
  tm_raster()+
tm_shape(squares_buffer)+
  tm_polygons()
#overlap fine

#extract land use for the kvatradnr
squares_buffer_landUse <- raster::extract(raster::raster(r),as_Spatial(squares_buffer),
                                          df=T,weights=TRUE)
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
squares_buffer <- bind_cols(squares_buffer,landuse_summary)

#plot to check
tm_shape(squares_buffer)+
  tm_fill("urban")          

tm_shape(squares_buffer)+
  tm_fill("agri_int")   

tm_shape(squares_buffer)+
  tm_fill("forest")  

saveRDS(squares_buffer,file="environ-data/squares_buffer_1km.rds")
saveRDS(squares_buffer,file="environ-data/squares_buffer_5km.rds")

### TRANSECTS ####

#read in dof squares
lines <- st_read(dsn = "data/Insects_and_TTT",
                   layer = "transects utm32")
plot(lines)#lines

#get land use data
ldir <- "C:/Users/db40fysa/Dropbox/DOF/landuseData/basemap03_2018"
#aggregated layer with all classes
r <- rast(paste(ldir,"lu_agg_2018.tif",sep="/"))
r <- aggregate(r, fact=2,fun="modal")
r[r==999999] <- NA
#sort(unique(values(r)))#corresponds with table 3.6!!
#table(values(r))
plot(r)
crs(r) <- "epsg:32632"

#overlay land map with lines
tm_shape(raster::raster(r))+
  tm_raster()+
  tm_shape(lines)+
  tm_lines()
#overlap fine

#small buffer on the lines
lines_buffer <- st_buffer(squares,dist=5)

#extract land use for the kvatradnr
lines_buffer_landUse <- raster::extract(raster::raster(r),as_Spatial(lines_buffer), 
                                        df=T, weights=TRUE)
names(lines_buffer_landUse)[2] <- "LandUse"
lines_buffer_landUse <- subset(lines_buffer_landUse, !is.na(LandUse))

#summarize per ID
lines_landUse_summary <- lines_buffer_landUse %>%
  dplyr::group_by(ID,LandUse) %>%
  dplyr::summarise(Length = sum(weight))

#see Table 3.6 in TR159.pdf for explanation of land use codes
table(lines_landUse_summary$LandUse)

#now get proportion of urban areas and agricultural areas
landuse_summary <- lines_landUse_summary %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(path = sum(Length[LandUse %in% 142000]),
                   road = sum(Length[LandUse %in% c(141000,150000)]), # or road or rail
                   forest = sum(Length[LandUse %in% 311000:312000]),
                   mapped = sum(Length[!LandUse %in% 800000]),
                   total = sum(Length))  

head(landuse_summary)


#combine with land use data
lines <- bind_cols(lines_buffer,landuse_summary) %>%
          as_tibble()

#plot to check
tm_shape(lines)+
  tm_lines("urban")          

tm_shape(lines)+
  tm_lines("agri_int")   

tm_shape(lines)+
  tm_lines("forest")  

saveRDS(lines,file="environ-data/lines.rds")
saveRDS(lines,file="environ-data/lines_5m.rds")
saveRDS(lines,file="environ-data/lines_20m.rds")
