### libraries ####

library(AHMbook)
library(tidyverse)
library(lubridate)
library(unmarked)
library(ggthemes)
library(cowplot)

### DOF data #####

source('processing.R', echo=TRUE)

### function ####

fitModel <- function(myspecies, output="ESW"){

#and just one species
dataS <- data %>% 
          filter(Species==myspecies) %>%
          select(-id)

#and use the right season for the species
dataS <- dataS %>%
          filter(type==bestSeason)

#join data - all years
allData <- left_join(info, dataS)

#add on environmental data
allData <- inner_join(allData, environData)

#get number of birds seen at each distance
distanceData <- allData %>% 
            dplyr::select(Year, X.0,X.1,X.2) %>%
            pivot_longer(starts_with("X."),names_to = "distance", values_to = "nu") %>%
            dplyr::filter(!is.na(nu)) %>%
            dplyr::mutate(distance = recode(distance, X.0 = "25", X.1 = "50", X.2 = "100")) %>%
            dplyr::group_by(Year, distance) %>%
            dplyr::summarise(total = sum(nu)) %>%
            dplyr::mutate(distance = factor(distance, levels = c("25", "50", "100"))) %>%
            ungroup() %>%
            add_column(Species = myspecies)
            
#organise data for the distance model

#was a quadrat number visited more than once? no
#summaryQuad <- allData %>%
#                    dplyr::group_by(kvadratnr) %>%
#                    dplyr::summarise(nu = length(unique(dato)))#all only visited once
#table(summaryQuad$nu)

#get number of individuals detected per site
# procData <- allData %>% 
#                 dplyr::select(kvadratnr,X.0,X.1,X.2) %>%
#                 pivot_longer(starts_with("X."),names_to="distance", values_to="nu") %>%
#                 dplyr::mutate(distance = 
#                                   recode(distance, X.0 = "25", X.1 = "50", X.2 = "100")) %>%
#                 dplyr::mutate(nu = ifelse(is.na(nu),0,nu),
#                               distance = as.numeric(as.character(distance)))
    

# #one per per detection
# procDataR <- procData %>%
#               filter(nu > 0) %>%
#               group_by(kvadratnr) %>%
#               slice(rep(1:n(), nu)) %>%
#               dplyr::select(kvadratnr, distance) %>%
#               ungroup()
# table(procDataR$distance)

### distance model ####

# library(Distance)
# ds_hn <- ds(procData, transect="line",key="hn")
# ds_hr <- ds(procData, transect="line",key="hr")
# ds_uni <- ds(procData, transect="line",key="unif")
# 
# plot(ds_hn)
# gof_ds(ds_hn)
# summarize_ds_models(ds_hn, ds_hr, ds_uni, output="plain")

### format data for unmarked ####

allData$line_pathroad
covariates <- data.frame(scale(allData[,c("skydaekke","regn","vind","Year",
                                          "lines_path","lines_road","lines_forest",
                                          "squares_forest","squares_agri_int",
                                          "squares_urban", "squares_agri_ext",
                                          "squares_freshwater")]))
covariates$fYear <- as.factor(allData$Year)

#distance data
temp <- allData[,c("X.0","X.1","X.2")]
temp[is.na(temp)] <- 0
colSums(temp)

### unmarked ####

#format data frame
unmarkDF <- unmarkedFrameDS(y=as.matrix(temp),
                            siteCovs= covariates,
                            dist.breaks=c(0,25,50,100), 
                            tlength=rep(1000,nrow(allData)),
                            unitsIn="m", 
                            survey="line")

#with covariates
#+ lines_forest
(fm1 <- distsamp(~ lines_path +lines_road 
                 ~ fYear + squares_forest + squares_agri_int + squares_urban + squares_agri_ext + squares_freshwater, 
                 data = unmarkDF, 
                 keyfun = "halfnorm",
                 output = "density",
                 unitsOut = "kmsq"))


#decide on output
if(output=="state"){
  
stateDF <- data.frame(Species = myspecies,
                       param = names(coef(fm1,type='state')),
                       coef = as.numeric(coef(fm1,type='state')),
                       coef_se = as.numeric(SE(fm1,type='state')))[-1,]
return(stateDF)

}else if (output=="detection"){
  
detectionDF <- data.frame(Species = myspecies,
                       param = names(coef(fm1,type='det')),
                       coef = as.numeric(coef(fm1,type='det')),
                       coef_se = as.numeric(SE(fm1,type='det')))[-1,]

return(detectionDF)

#use model to predict sigma
}else if(output=="ESW"){

log_sigma <- predict(fm1, type="det", newdata = data.frame(lines_path = min(covariates$lines_path),
                                              lines_road = min(covariates$lines_road),
                                              lines_forest = covariates$lines_forest))
#get effective strip width
log_sigma$ESW <- sapply(log_sigma$Predicted, function(x){
  integrate(gxhn, 0, 101, sigma = x)$value
})
#head(log_sigma)#sigma is 38.2, ESW is 47 which makes sense
#backTransform(fm1, type="det")  #same as what is in predict             
#sqrt((pi * 38.51699^2)/2)

#add on species
log_sigma$Species <- myspecies
return(log_sigma)

}else if(output=="simple_ESW"){
  
  #null model
  (fm0 <- distsamp(~ 1
                   ~1, 
                   data = unmarkDF, 
                   keyfun = "halfnorm",
                   output = "density",
                   unitsOut = "kmsq"))
  
  log_sigma <- predict(fm0, type="det",
                       newdata = data.frame(lines_path = min(covariates$lines_path),
                                            lines_road = min(covariates$lines_road),
                                            lines_forest = covariates$lines_forest))
  
  #get effective strip width
  log_sigma$ESW <- sapply(log_sigma$Predicted, function(x){
    integrate(gxhn, 0, 101, sigma = x)$value
  })
  #head(log_sigma)#sigma is 38.2, ESW is 47 which makes sense
  #backTransform(fm1, type="det")  #same as what is in predict             
  #sqrt((pi * 38.51699^2)/2)
  
  #add on species
  log_sigma$Species <- myspecies
  return(log_sigma[1,])
  
} else if(output=="density"){
  
  newdata <- covariates
  newdata$lines_path <- min(newdata$lines_path)
  newdata$lines_road <- min(newdata$lines_road)
  
  densityDF <- predict(fm1, type="state", newdata = newdata)
  
  #add on species
  densityDF$Species <- myspecies
  return(densityDF)

} else if(output=="distanceData"){
  
  return(distanceData)
  
}

}


### plot distance data ####

distanceData <- lapply(species, function(x){
  fitModel(x,output="distanceData")}) %>%
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

### ESWs #####

#simple ESW
simpleDF <- lapply(species, function(x){
  fitModel(x,output="simple_ESW")}) %>%
  reduce(rbind)

#predicted ESW including covariates
covDF <-lapply(species, function(x){
  fitModel(x,output="ESW")}) %>%
  reduce(rbind)
  
#combined and add on mean ESW
outputDF <- covDF %>%
                left_join(.,simpleDF, by="Species", suffix=c("_C","_N")) %>%
                group_by(Species) %>%
                mutate(meanESW = mean(ESW_C)) %>%
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
              summarise(ESW = median(ESW_C)) %>%
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
            Diet.5Cat + flockSize + Forest_sp,
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

