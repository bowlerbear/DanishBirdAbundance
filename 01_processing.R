### libraries ####

library(tidyverse)
library(lubridate)
library(ggthemes)
library(cowplot)
library(sf)

### DOF data #####

#read in bird observation data - counts of birds using distance sampling
data <- read.csv("data/TTT/ttt_data.csv",sep=";")
data$total <- data$X.0 + data$X.1 + data$X.2
data$Date <- as.Date(data$dato)
data$Year <- year(data$Date)

#read in transect info
info <- read.csv("data/TTT/ttt_info.csv",sep=";")
info$Date <- as.Date(info$dato)
info$Year <- year(info$Date)

### subset data ####

#start with just the last year of data??
#infoS <- info %>%
#  filter(Year==2016) %>%
#  select(-id)

#not winter
data <- data %>%
          dplyr::filter(type!="vinter") 

### sort species list ####

#get list from the Danish guys
danishSpecies <- read.csv("~/Dropbox/DOF/info/species_max_season_comments.csv")

#sort (sub)species to be combined, as suggested by the DANES
combineSpecies <-  danishSpecies %>%
                    filter(grepl(";",EURING_include))

toCombine <- lapply(1:nrow(combineSpecies), function(x){
  
  allEURINGS <- strsplit(as.character(combineSpecies$EURING_include[x]),";")[[1]]
  Species <- unique(danishSpecies$Scientific_name[danishSpecies$EURING_subspecies %in% allEURINGS])
  finalSpecies <- combineSpecies$Scientific_name[x]
  EURING <- combineSpecies$EURING_present_as[x]
  data.frame(finalSpecies,EURING,Species)
  
}) %>% do.call(rbind,.)

#add this info on combining names to the original datasets
data$finalSpecies <- toCombine$finalSpecies[match(data$latin,toCombine$Species)]
data$Species <- data$latin  
data$Species[!is.na(data$finalSpecies)] <- data$finalSpecies[!is.na(data$finalSpecies)]

#add info on whether to estimate that species
data$use <- danishSpecies$Produce_estimate[match(data$Species, danishSpecies$Scientific_name)]
data <- data %>% filter(use=="yes")
sort(unique(data$Species))

#add info on best season to use for each species
data$bestSeason <- danishSpecies$Stage_use[match(data$Species, danishSpecies$Scientific_name)]
table(data$bestSeason)

### select species ####

#when was each species most often seen?
speciesSeason <- data %>%
                    dplyr::group_by(Species,type) %>%
                    dplyr::summarise(total=sum(total)) %>%
                    dplyr::group_by(Species) %>%
                    dplyr::summarise(type=type[total==max(total)], 
                                     total = sum(total)) %>%
                    dplyr::arrange(desc(total))

table(speciesSeason$type)
#in the end, use the expert-decided bestSeason above
                      
#select species to fill model to
species <- c("Alauda arvensis", "Perdix perdix", 
               "Passer domesticus", "Passer montanus",
               "Fringilla coelebs","Hirundo rustica")
myspecies <- species[5]

#or most common species
species <- speciesSeason$Species[speciesSeason$total>500]
length(species)

#or rarer species
species <- speciesSeason$Species[speciesSeason$total<500 & speciesSeason$total>50]
length(species)

#or passerines
traits <- readRDS("traits/traits.rds")
species <- speciesSeason$Species[speciesSeason$total>50 & 
                                   speciesSeason$Species %in% traits$Species[traits$PassNonPass=="Passeriformes"]]

### flock size data ###
flockSize <- data %>%
              group_by(Species) %>%
              summarize(flockSize = median(total[total!=0]))

### environdata ####

#line data - factor affecting detectability
lines <- readRDS("environ-data/lines_5m.rds")
names(lines)[4:8] <- sapply(names(lines)[4:8], function(x){paste0("lines_",x)})
lines$lines_pathroad <- lines$lines_path + lines$lines_road
lines$lines_pathroad <- log(lines$lines_pathroad/lines$lines_mapped+0.01)
lines$lines_path <- log(lines$lines_path/lines$lines_mapped+0.01)
lines$lines_road <- log(lines$lines_road/lines$lines_mapped+0.01)
lines$lines_forest <- log(lines$lines_forest/lines$lines_mapped+0.01)

#buffer data - factors affecting abundance
squares <- readRDS("environ-data/squares_buffer_1km.rds")
names(squares)[4:13] <- sapply(names(squares)[4:13], function(x){paste0("squares_",x)})
squares$squares_forest <- log(squares$squares_forest/squares$squares_mapped+0.01)
squares$squares_agri_int <- squares$squares_agri_int/squares$squares_mapped
squares$squares_urban <- log(squares$squares_urban/squares$squares_mapped+0.01)
squares$squares_agri_ext <- log(squares$squares_agri_ext/squares$squares_mapped+0.01)
squares$squares_freshwater <- log(squares$squares_freshwater/squares$squares_mapped+0.01)

environData <- inner_join(lines[,c("kvadratnr","lines_path","lines_road","lines_pathroad",
                                   "lines_forest")],
                      squares[,c("kvadratnr","squares_forest","squares_agri_int",
                                 "squares_urban", "squares_agri_ext",
                                 "squares_freshwater")])

#pairs(environData)

#compare density estimates with and without accounting for covariates
#cuckoo maybe only seen
#fw data
#people told to focus on the strip
#smaller species harder to distinguish when further away

#model variation in ESW due to traits:
#size
#flocking behaviour
#sound clues
#foraging stratum

### exact transect lengths ######

# transects32 <- st_read(dsn = "data/TTT",
#                      layer = "transects utm32")
# 
# transects33 <- st_read(dsn = "data/TTT",
#                        layer = "transects utm32") %>%
#                 st_transform(.,crs = st_crs(transects32))
# all(transects32$kvadratnr==transects33$kvadratnr)
# #they are the same...
# 
# #get exact lengths
# transects33$length <- st_length(transects33)
# summary(transects33$length)
# #all pretty much 1km
