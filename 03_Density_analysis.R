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

#basic observed stats

speciesSummary <- data %>%
  filter(type==bestSeason) %>%
  group_by(Species) %>%
  summarise(obs_count = sum(total))

### choose models #####

myModels <- readRDS("outputs/models_passerines_lines_path_road.rds")

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

  #cap the maximum number observation at the 99% quantile
  top99 <- as.numeric(round(quantile(allData$total,0.99)))
  allData$total <- ifelse(allData$total > top99,
                          top99, allData$total)
  hist(allData$total)
  
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
output1 <- species %>%
  map_dfr(.,fitGAM)
saveRDS(output,file="output.rds")

#### upscaling ###
output1 <- readRDS("output.rds")

#add on fraction surveyed and ESW
output1$surveyFraction <- nullDF$surveyFraction[match(output1$Species, nullDF$Species)]
output1$ESW <- nullDF$ESW[match(output1$Species, nullDF$Species)]

#predict total number per square
output1$scaled_preds <- output1$preds/output1$surveyFraction
output1$scaled_preds_lower <- (output1$preds - 2*output1$preds_sd)/output1$surveyFraction
output1$scaled_preds_upper <- (output1$preds + 2*output1$preds_sd)/output1$surveyFraction

#add on observed count
output1 <- inner_join(output1,speciesSummary, by="Species")

#total per species
(outputSummary <- output1 %>%
  group_by(Species) %>%
  summarise(pred_count =sum(preds),
            total_count = sum(scaled_preds),
            total_count_lower = sum(scaled_preds_lower),
            total_count_upper = sum(scaled_preds_upper),
            obs_count = mean(obs_count),
            ESW = mean(ESW),
            surveyFraction = mean(surveyFraction)) %>%
  arrange(desc(total_count))) 

#remove outliders
outputSummary <- outputSummary %>%
                  filter(!Species %in% c("Turdus pilaris",
                                         "Loxia curvirostra"))
#plotting
ggplot(outputSummary) + 
  geom_text(aes(x=pred_count,y=total_count,label=Species))

ggplot(outputSummary) + 
  geom_text(aes(x=obs_count,y=total_count,label=Species),size=rel(2))

outputSummary %>%
  filter(obs_count<4000) %>%
  ggplot() + 
  geom_text(aes(x=obs_count,y=total_count,label=Species),size=rel(2))
  
ggplot(outputSummary)+
  geom_col(aes(x = fct_reorder(Species, total_count), y = total_count, 
               fill = ESW)) +
  xlab("Species") + ylab("Estimated total population size") +
  coord_flip() +
  theme_bw()

### uncertainty ##########

#bootstrapping to estimate uncertainty?
#delta method?
#https://examples.distancesampling.org/Distance-variance/variance-distill.html
#http://workshops.distancesampling.org/standrews-2019/intro/lectures/BlockD-precision-poststrat.pdf
#Horvitz-Thompson-like estimator

### irregular squares ####

#### upscaling ####