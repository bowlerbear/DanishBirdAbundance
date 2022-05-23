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

#based in environmenta data in the 1 x 1 km square
saveRDS(myModels, file="outputs/models_passerines_linespathroad.rds")#no year term
saveRDS(myModels, file="outputs/models_passerines_lines_path_road.rds")#sep terms

#based on environmental data in the 1 x 200 m transect area


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

source('00_functions.R', echo=TRUE)
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


### jacknifing ####################################

library(broom)
library(tidyverse)

# do a leave-one-out type thing
# where each species is left out at a time and
# then their ESW is predicted
leave_one_out_function <- function(species_name){
  
  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)
  
  # I removed "habitat from the model because it is tricky because
  # there are two levels where only one species is that habitat
  mod <- lm(ESW ~ log(Mass) + ForStrat.ground + 
              Trophic.Level + Primary.Lifestyle + Trophic.Niche,
            data=dat_filtered)
  
  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)
  
  predicted_ESW <- predict(mod, new_dat)
  
  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)
  
  return(loo_summary)
  
}

leave_one_out_results <- bind_rows(lapply(unique(traitsDF$Species), leave_one_out_function))

ggplot(leave_one_out_results, aes(x=observed_ESW, y=predicted_ESW))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Observed ESW")+
  ylab("Predicted ESW")+
  geom_smooth(method="lm")

summary(lm(observed_ESW ~ predicted_ESW, data=leave_one_out_results))

# but how many species can we leave out at a time and still get a good predictive accuracy of the model?
random_sampling_function <- function(sample_size){
  
  boot_fun <- function(draw_number){
    
    dat_filtered <- traitsDF %>%
      sample_n(nrow(traitsDF)-sample_size)
    
    mod <- lm(ESW ~ log(Mass) + ForStrat.ground +  
                Trophic.Level + Primary.Lifestyle + Trophic.Niche,
              data=dat_filtered)
    
    new_dat <- traitsDF %>%
      dplyr::filter(!Species %in% dat_filtered$Species)
    
    predicted_ESW <- predict(mod, new_dat)
    
    loo_summary <- data.frame(Species=new_dat$Species,
                              observed_ESW=new_dat$ESW,
                              predicted_ESW=predicted_ESW) %>%
      mutate(draw=draw_number)
    
    return(loo_summary)
    
  }
  
  boot_results <- bind_rows(lapply(c(1:50), boot_fun)) %>%
    mutate(sample=sample_size) %>%
    mutate(percent_missing=(sample/nrow(traitsDF))*100)
 
  
}

lapply_with_error <- function(X,FUN,...){    
  lapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                      error=function(e) NULL))
}

# adding the lapply with error
# allows to 'skip'
# the possibility when something is left out
# that has a level not in the model
# e.g., primary lifestyle this happens
random_sampling_results <- bind_rows(lapply_with_error(c(1:60), random_sampling_function))

# quick crack at a potential visualization
# not sure what I'm doing
ggplot(random_sampling_results, aes(x=percent_missing, y=predicted_ESW))+
  geom_point()

# that didn't look too good...

# try again?
random_sampling_results %>%
  mutate(percent_missing=as.numeric(percent_missing)) %>%
  nest(data = -percent_missing) %>% 
  mutate(
    fit = map(data, ~ lm(predicted_ESW ~ observed_ESW, data = .x)),
    tidied = map(fit, tidy),
    glanced=map(fit, glance)
  ) %>% 
  unnest(glanced) %>%
  ggplot(., aes(x=percent_missing, y=r.squared))+
  geom_jitter()


# so that kind of makes sense
# but with the lapply with error....
# by the time you have the majority of species being 'left-out'
# then the levels are gone so it isn't working.
# just as a test case I'm gonna copy and paste the above
# but just use body mass (as we know we have this for all species and don't have to worry about levels)
# but how many species can we leave out at a time and still get a good predictive accuracy of the model?
random_sampling_function_body <- function(sample_size){
  
  boot_fun <- function(draw_number){
    
    dat_filtered <- traitsDF %>%
      sample_n(nrow(traitsDF)-sample_size)
    
    mod <- lm(ESW ~ log(Mass),
              data=dat_filtered)
    
    new_dat <- traitsDF %>%
      dplyr::filter(!Species %in% dat_filtered$Species)
    
    predicted_ESW <- predict(mod, new_dat)
    
    loo_summary <- data.frame(Species=new_dat$Species,
                              observed_ESW=new_dat$ESW,
                              predicted_ESW=predicted_ESW) %>%
      mutate(draw=draw_number)
    
    return(loo_summary)
    
  }
  
  boot_results <- bind_rows(lapply(c(1:50), boot_fun)) %>%
    mutate(sample=sample_size) %>%
    mutate(percent_missing=(sample/nrow(traitsDF))*100)
  
  
}

# adding the lapply with error
# allows to 'skip'
# the possibility when something is left out
# that has a level not in the model
# e.g., primary lifestyle this happens
random_sampling_results <- bind_rows(lapply_with_error(c(1:60), random_sampling_function_body))

# try again?
random_sampling_results %>%
  mutate(percent_missing=as.numeric(percent_missing)) %>%
  nest(data = -percent_missing) %>% 
  mutate(
    fit = map(data, ~ lm(predicted_ESW ~ observed_ESW, data = .x)),
    tidied = map(fit, tidy),
    glanced=map(fit, glance)
  ) %>% 
  unnest(glanced) %>%
  ggplot(., aes(x=percent_missing, y=r.squared))+
  geom_jitter()

# Ah, that's more what I was expecting!
# Kinda cool... I think?


### det covariates ####

#dextract coefficients for the etection covariates in the model
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