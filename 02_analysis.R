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

### plot distance data ####

distanceData <- lapply(species, function(x){
  getDistanceData(x)
  }) %>%
  reduce(rbind)

ggplot(distanceData) + 
   geom_col(aes(x = distance, y = total, fill = Year), position="stack")+
   facet_wrap(~Species,scales="free")+
   scale_fill_viridis_c()

distanceData$Common <- data$english[match(distanceData$Species, data$Species)]
sort(unique(distanceData$Common))

saveRDS(distanceData, file="outputs/distanceData_passerines.rds")

### problem species ####

probs_set1 <- c("Anser anser", "Aythya fuligula", "Buteo buteo", "Chroicocephalus ridibundus", 
                "Fulica atra", "Larus argentatus", "Larus canus", "Anas platyrhynchos", 
                "Phasianus colchicus", "Vanellus vanellus")

probs_set2 <- c("Ardea cinerea", "Circus aeruginosus", "Columba livia", 
                "Columba oenas", "Corvus corax", "Corvus corone", "Cuculus canorus", 
                "Cygnus olor", "Falco tinnunculus", "Gallinula chloropus","Larus fuscus", 
                "Luscinia luscinia", "Phalacrocorax carbo", "Podiceps cristatus", 
                "Tadorna tadorna", "Tringa totanus")

select_species <- species[!species %in% probs_set1]
select_species <- species[!select_species %in% probs_set2]

### fit models ###

myModels <- lapply(species, function(x){
  fitModel(x)
})

saveRDS(myModels, file="outputs/models_passerines_fYear_linespathroad.rds")

### extract ####

### ESWs #####

#simple ESW
simpleDF <- lapply(myModels, function(x){
  extractModel(x,"fm0",output="ESW")
  }) %>%
  reduce(rbind)

#predicted ESW including covariates
covDF <-lapply(species, function(x){
  extractModel(x,"fm1",output="ESW")
  }) %>%
  reduce(rbind)
  
#combined and add on mean ESW
outputDF <- covDF %>%
                left_join(.,simpleDF, by="Species", suffix=c("_C","_N")) %>%
                group_by(Species) %>%
                mutate(meanESW = median(ESW_C)) %>%
                ungroup() %>%
                mutate(Species = fct_reorder(factor(Species), meanESW))
saveRDS(outputDF, file="outputs/outputDF_passerines.rds")


#plotting the ESW of each species
ggplot(outputDF)+
  geom_boxplot(aes(x=Species, y = ESW_C))+
  geom_point(aes(x=Species, y = ESW_N),colour="red")+
  coord_flip()+
  theme_few()+
  ylab("Effective strip width (m)")

### traits ####

#effects of traits
outputDF <- readRDS("outputs/outputDF_passerines.rds")
traits <- readRDS("traits.rds")
traitsDF <- outputDF %>%
              group_by(Species) %>%
              summarise(ESW = median(ESW_N)) %>%
              inner_join(., traits, by="Species") %>% 
              inner_join(.,flockSize, by="Species")

theme_set(theme_bw())
g1 <- qplot(log(BodyMass.Value), ESW, data=traitsDF) + stat_smooth(method="lm")
g2 <- qplot(Diet.5Cat, ESW, data=traitsDF, geom="boxplot")
qplot(Diet.5Cat, ESW, data=traitsDF)
g3 <- qplot(ForStrat.ground, ESW, data=traitsDF) + stat_smooth(method="lm")
g4 <- qplot(ForStrat.aerial, ESW, data=traitsDF) + stat_smooth(method="lm")
qplot(flockSize, ESW, data=traitsDF)
qplot(Forest_sp, ESW, data=traitsDF)

plot_grid(g1,g2,g3,g4,ncol=2)

qplot(Habitat, ESW, data=traitsDF, geom="boxplot")
qplot(factor(Habitat.Density), ESW, data=traitsDF, geom="boxplot")
qplot(Trophic.Level, ESW, data=traitsDF, geom="boxplot")
qplot(Trophic.Niche, ESW, data=traitsDF, geom="boxplot")
qplot(Primary.Lifestyle, ESW, data=traitsDF, geom="boxplot")

#models
lm1 <- lm(ESW ~ ForStrat.ground + ForStrat.aerial + log(BodyMass.Value) + 
            Diet.5Cat + flockSize,
          data=traitsDF)

lm1 <- lm(ESW ~ log(Mass) + Trophic.Level + Primary.Lifestyle + Habitat.Density,
          data=traitsDF)

summary(lm1)

### det covariates ####

#trouble with: 
#"Passer montanus"     "Riparia riparia"     "Aegithalos caudatus" "Certhia familiaris" 

#detection
detectionDF <-lapply(species[-c(4,46,55,56)], function(x){
  fitModel(x,output="detection")}) %>%
  reduce(rbind)

saveRDS(detectionDF, file="outputs/detectionDF_passerines.rds")

detectionDF %>%
  filter(!Species %in% c("Spinus spinus","Phoenicurus ochruros",
                         "Corvus corone")) %>%
  ggplot()+
  geom_crossbar(aes(x=Species,y=coef,ymin=coef-coef_se,ymax=coef+coef_se))+
  facet_wrap(~param)+coord_flip()+theme_few()+
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Effect on ESW")


### state covariates ####
#state
stateDF <-lapply(species, function(x){
  fitModel(x,output="state")}) %>%
  reduce(rbind)

stateDF %>%
  #filter(!Species %in% c("Riparia riparia","Columba livia","Cuculus canorus")) %>%
  ggplot()+
  geom_crossbar(aes(x=Species,y=coef,ymin=coef-coef_se,ymax=coef+coef_se))+
  facet_wrap(~param)+coord_flip()+theme_few()+
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Effect on abundance")


### density predictions ####

densityDF <-lapply(species[-5], function(x){
  fitModel(x,output="density")}) %>%
  reduce(rbind)

densityDF %>%
  group_by(Species) %>%
  mutate(medPred = median(Predicted)) %>%
  ungroup() %>%
  mutate(Species = fct_reorder(factor(Species), medPred)) %>%
  ggplot()+
  geom_violin(aes(x=Species,y=Predicted), fill="blue", alpha=0.5)+
  coord_flip()+
  scale_y_log10()+
  ylab("Predicted density km2")+
  theme_few()

### upscaling ####
library(sf)

#area of squares
squares32 <- st_read(dsn = "data/Insects_and_TTT", layer = "transect squares utm32")
squares33 <- st_read(dsn = "data/Insects_and_TTT", layer = "transect squares utm33") %>%
              st_transform(.,st_crs(squares32))
squares <- bind_rows(squares32, squares33)
squares$area <- st_area(squares)
sum(squares$area)

#area of Denmark
Denmark <- raster::getData('GADM', country='DNK', level=0) %>%
            st_as_sf() %>%
            st_transform(., st_crs(squares32)) %>%
            add_column(area = st_area(.))

sum(Denmark$area)

#fraction surveyed
sum(squares$area)/sum(Denmark$area)
#29%

