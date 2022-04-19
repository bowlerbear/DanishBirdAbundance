### libraries ####

library(AHMbook)
library(tidyverse)
library(lubridate)
library(unmarked)
library(ggthemes)
library(cowplot)

### DOF data #####

source('00_functions.R', echo=TRUE)
source('01_processing.R', echo=TRUE)

### choose models #####

myModels <- readRDS("outputs/models_passerines_linespathroad.rds")

### ESWs #####

#fm0- null model
nullDF <- lapply(myModels, function(x){
  extractModel(x,"fm0", output="ESW")
}) %>%
  reduce(rbind) %>%
  add_column(type="null")

#fm3 - no detection covariates
simpleDF <- lapply(myModels, function(x){
  extractModel(x,"fm3", output="ESW")
}) %>%
  reduce(rbind) %>%
  add_column(type="no covariates")

#### fraction surveyed ####

#fraction of square effectively surveyed

nullDF$transectArea <- 1000 * nullDF$ESW * 2
squareArea <- 1000 * 1000
nullDF$surveyFraction <- nullDF$transectArea/squareArea
summary(nullDF$surveyFraction)

#### total sampled area ####

library(sf)

#area of squares
squares32 <- st_read(dsn = "data/TTT", layer = "transect squares utm32")
squares33 <- st_read(dsn = "data/TTT", layer = "transect squares utm33") %>%
              st_transform(.,st_crs(squares32))
squares <- bind_rows(squares32, squares33) %>%
              filter(kvadratnr %in% data$kvadratnr) %>%
              filter(!duplicated(kvadratnr))
              
squares$area <- st_area(squares)
sum(squares$area)
nrow(squares)

#area of Denmark
Denmark <- raster::getData('GADM', country='DNK', level=0) %>%
            st_as_sf() %>%
            st_transform(., st_crs(squares32)) %>%
            add_column(area = st_area(.))

sum(Denmark$area)

#fraction surveyed
(sampledRatio <- sum(squares$area)/sum(Denmark$area))
#0.04%

### 1km match squares - gam ####

#grid for whole Denmark
gridData <- readRDS("environ-data/grid_environdata.rds")
gridData$squares_forest <- log(gridData$forest/gridData$mapped+0.01)
gridData$squares_agri_int <- gridData$agri_int/gridData$mapped
gridData$squares_urban <- log(gridData$urban/gridData$mapped+0.01)
gridData$squares_agri_ext <- log(gridData$agri_ext/gridData$mapped+0.01)
gridData$squares_freshwater <- log(gridData$freshwater/gridData$mapped+0.01)

fitGAM <- function(myspecies, plot=FALSE){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    select(-id)
  
  #and use the right season for the species
  dataS <- dataS %>%
    mutate(type = as.character(type)) %>%
    mutate(bestSeason=as.character(bestSeason)) %>%
    dplyr::filter(type==bestSeason)
  
  #join data - all years
  allData <- left_join(dataS,info)
  
  #add on environmental data
  allData <- inner_join(allData, environData)
  
  #add coordinates
  allData <- addCoords(allData)

  #fit gam  
  library(brms)
  gam1 <- brm(bf(total ~ skydaekke + regn + vind + 
                squares_forest + squares_agri_int + squares_urban +
                squares_agri_ext + squares_freshwater  + s(X,Y)), 
              family=poisson,
              data=allData)
  
  #predict to whole extent
  temp <- predict(gam1, newdata = expand_grid(gridData,
                                                       skydaekke = 1,
                                                       regn = 1,
                                                       vind = 1))
  gridData$preds <- temp[,"Estimate"]
  gridData$preds_sd <- temp[,"Est.Error"]
  gridData$Species <- myspecies
  
  if(plot==TRUE){
  g1 <- ggplot(gridData)+
    geom_tile(aes(x=X,y=Y, fill=preds))+
    scale_fill_viridis_c(trans="log10")
    
  g2 <- ggplot(allData)+
    geom_point(aes(x=X,y=Y, colour=total))+
    scale_colour_viridis_c(trans="log10")
  
  cowplot::plot_grid(g1,g2)
  }
  
  return(gridData)
  
  }

#get predictions for each  species
output1 <- species[1:10] %>%
  map_dfr(.,fitGAM)
saveRDS(output1,file="output1.rds")


#### upscaling ####


#add on fraction surveyed


### irregular squares ####


#### upscaling ####


