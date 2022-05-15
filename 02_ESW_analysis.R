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

### distance data ####

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

#these species do not show the classic distance decay

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

saveRDS(myModels, file="outputs/models_passerines_linespathroad.rds")#no year term
saveRDS(myModels, file="outputs/models_passerines_lines_path_road.rds")#sep terms

### pick model ####

#fm0 = null, fm1 = cov, fm2 = null state, fm3 = null detection

myModels <- readRDS("outputs/models_passerines_linespathroad.rds")
myModels <- readRDS("outputs/models_passerines_lines_path_road.rds")

### ESWs #####

#fm0- null model
nullDF <- lapply(myModels, function(x){
  extractModel(x,"fm0", output="ESW")
}) %>%
  reduce(rbind) %>%
  add_column(type="null")

summary(nullDF$ESW)

#fm3 - no detection covariates
simpleDF <- lapply(myModels, function(x){
  extractModel(x,"fm3", output="ESW")
}) %>%
  reduce(rbind) %>%
  add_column(type="no covariates")

summary(simpleDF$ESW)

# with detection covariates - single terms
combDF <- lapply(myModels, function(x){
  extractModel(x,"fm1", output="ESW", sep=FALSE)
}) %>%
  reduce(rbind) %>%
  add_column(type="combined terms")

# with detection covariates - separate terms
myModels <- readRDS("outputs/models_passerines_lines_path_road.rds")
sepDF <- lapply(myModels, function(x){
  extractModel(x,"fm1", output="ESW")
}) %>%
  reduce(rbind) %>%
  add_column(type="separate terms")

#### compare different model ESW #####

simpleDF %>%
  inner_join(.,nullDF, by="Species") %>%
  ggplot() +
  geom_point(aes(x = ESW.x, y = ESW.y)) +
  geom_abline()
#almost the same

simpleDF %>%
    inner_join(.,sepDF, by="Species") %>%
    ggplot() +
    geom_point(aes(x = ESW.x, y = ESW.y)) +
    geom_abline()

nullDF %>%
  inner_join(.,combDF, by="Species") %>%
  ggplot() +
  geom_point(aes(x = ESW.x, y = ESW.y)) +
  geom_abline()

sepDF %>%
  inner_join(.,combDF, by="Species") %>%
  ggplot() +
  geom_point(aes(x = ESW.x, y = ESW.y)) +
  geom_abline()

### traits ####

myModels <- readRDS("outputs/models_passerines_linespathroad.rds")

simpleDF <- lapply(myModels, function(x){
  extractModel(x,"fm3", output="ESW")
}) %>%
  reduce(rbind) %>%
  add_column(type="no covariates")

summary(simpleDF$ESW)

#effects of traits
traits <- readRDS("traits/traits.rds")
traitsDF <- nullDF %>%
              group_by(Species) %>%
              summarise(ESW = median(ESW)) %>%
              ungroup() %>%
              inner_join(., traits, by="Species") %>% 
              inner_join(.,flockSize, by="Species")

theme_set(theme_bw())
g1 <- qplot(log(BodyMass.Value), ESW, data=traitsDF) + stat_smooth(method="lm")
g2 <- qplot(Diet.5Cat, ESW, data=traitsDF, geom="boxplot")
qplot(Diet.5Cat, ESW, data=traitsDF)
g3 <- qplot(ForStrat.ground, ESW, data=traitsDF) + stat_smooth(method="lm")
qplot(flockSize, ESW, data=traitsDF)

g4 <- qplot(Habitat, ESW, data=traitsDF, geom="boxplot") + coord_flip()
g5 <- qplot(factor(Habitat.Density), ESW, data=traitsDF, geom="boxplot")
g6 <- qplot(Trophic.Level, ESW, data=traitsDF, geom="boxplot")
g7 <- qplot(Trophic.Niche, ESW, data=traitsDF, geom="boxplot")
g8 <- qplot(Primary.Lifestyle, ESW, data=traitsDF, geom="boxplot")

plot_grid(g1,g3,
          g7,g6,
          g4,g8,nrow=3)

#models
summary(lm(ESW ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level + Primary.Lifestyle + Trophic.Niche,
            data=traitsDF))
#75%

### det covariates ####

#detection covariates in the model
detectionDF <-lapply(myModels, function(x){
  extractModel(x,"fm1",output="detection")
}) %>%
  reduce(rbind)

length(detectionDF$coef)
sum(abs(detectionDF$coef)>=2)

detectionDF %>%
  filter(abs(coef) < 3) %>%
  ggplot()+
  geom_crossbar(aes(x=Species,y=coef,ymin=coef-coef_se,ymax=coef+coef_se))+
  facet_wrap(~param)+coord_flip()+theme_few()+
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Effect on ESW")

### using ds package ####

deer.df <- ds(sikadeer, key="hn", truncation="10%", convert_units = conversion.factor)
plot(deer.df, main="Half normal detection function")

### end ####