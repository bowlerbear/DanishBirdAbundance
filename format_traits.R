library(tidyverse)

### get species ####
species <- read.csv("~/Dropbox/DOF/info/species_max_season_comments.csv") %>%
            filter(Produce_estimate=="yes") %>%
            filter(str_detect(Scientific_name," ")) %>%
            pull("Scientific_name") %>%
            unique() %>% 
            sort()

#write.table(data.frame(Species=species), file="DOF_species.txt",
#            row.names=FALSE, sep="\t")

### elton traits ####
elton <- read.delim("~/Dropbox/DOF/traits/Elton/BirdFuncDat.txt",
                    as.is=TRUE)
#check species names

species[!species %in% elton$Scientific]
#[1] "Acanthis flammea/A. cabaret" "Chloris chloris"             "Chroicocephalus ridibundus"  "Coloeus monedula"           
#[5] "Corvus cornix"               "Curruca communis"            "Curruca curruca"             "Cyanistes caeruleus"        
#[9] "Emberiza calandra"           "Larus canus canus"           "Linaria cannabina"           "Lophophanes cristatus"      
#[13] "Mareca strepera"             "Periparus ater"              "Poecile montanus"            "Poecile palustris"          
#[17] "Saxicola rubicola"           "Spatula clypeata"            "Spatula querquedula"         "Spinus spinus"

elton$Scientific[elton$Scientific %in% c("Carduelis flammea")] <- "Acanthis flammea/A. cabaret" 
elton$Scientific[elton$Scientific %in% c("Carduelis chloris")] <- "Chloris chloris" 
elton$Scientific[elton$Scientific %in% c("Larus ridibundus")] <- "Chroicocephalus ridibundus" 
elton$Scientific[elton$Scientific %in% c("Corvus monedula")] <- "Coloeus monedula"
elton$Scientific[elton$Scientific %in% c("Sylvia communis")] <- "Curruca communis"
elton$Scientific[elton$Scientific %in% c("Sylvia curruca")] <- "Curruca curruca"
elton$Scientific[elton$Scientific %in% c("Parus caeruleus")] <- "Cyanistes caeruleus"
elton$Scientific[elton$Scientific %in% c("Miliaria calandra")] <- "Emberiza calandra"
elton$Scientific[elton$Scientific %in% c("Carduelis cannabina")] <- "Linaria cannabina"
elton$Scientific[elton$Scientific %in% c("Parus cristatus")] <- "Lophophanes cristatus"
elton$Scientific[elton$Scientific %in% c("Anas strepera")] <- "Mareca strepera"
elton$Scientific[elton$Scientific %in% c("Parus ater")] <- "Periparus ater"
elton$Scientific[elton$Scientific %in% c("Parus montanus")] <- "Poecile montanus"
elton$Scientific[elton$Scientific %in% c("Parus palustris")] <- "Poecile palustris"
elton$Scientific[elton$Scientific %in% c("Saxicola torquatus")] <- "Saxicola rubicola"
elton$Scientific[elton$Scientific %in% c("Anas clypeata")] <- "Spatula clypeata"
elton$Scientific[elton$Scientific %in% c("Anas querquedula")] <- "Spatula querquedula"
elton$Scientific[elton$Scientific %in% c("Carduelis spinus")] <- "Spinus spinus"

#add on old sub species for crows
elton2 <-  subset(elton, Scientific=="Corvus corone")
elton2$Scientific <- "Corvus cornix"
elton <- rbind(elton, elton2)
  
species[!species %in% elton$Scientific]

### coreys data ###

# flock size??

### join ####

myElton <- subset(elton, Scientific %in% species) %>%
            select(Scientific,BLFamilyLatin,Diet.5Cat, 
                   ForStrat.ground,ForStrat.aerial,
                   Nocturnal,BodyMass.Value,PassNonPass) %>%
            rename(Species = "Scientific")



#### habitat data #####
tdir <- "~/Dropbox/DOF/traits/Storchova/doi_10.5061_dryad.n6k3n__v1"
habitat <- read.delim(paste(tdir,"Life-history characteristics of European birds.txt",sep="/"),
                      as.is=T)

species[!species %in% habitat$Species]
habitat$Species[grepl("flammea",habitat$Species)]

#fix names
habitat$Species[habitat$Species %in% c("Acanthis flammea")] <- "Acanthis flammea/A. cabaret"
habitat$Species[habitat$Species %in% c("Larus ridibundus")] <- "Chroicocephalus ridibundus"
habitat$Species[habitat$Species %in% c("Corvus monedula")] <- "Coloeus monedula"
habitat$Species[habitat$Species %in% c("Sylvia communis")] <- "Curruca communis"
habitat$Species[habitat$Species %in% c("Sylvia curruca")] <- "Curruca curruca"
habitat$Species[habitat$Species %in% c("Saxicola torquatus")] <- "Saxicola rubicola"

habitat2 <-  subset(habitat, Species=="Corvus corone")
habitat2$Species <- "Corvus cornix"
habitat <- rbind(habitat, habitat2)

alltraits <- inner_join(myElton, habitat, by="Species")

### clean habitat ####

alltraits$Forest <- apply(alltraits[,c("Deciduous.forest","Coniferous.forest",
                                       "Woodland")],1,max)
alltraits$Open <- apply(alltraits[,c("Tundra", "Grassland",
                                     "Mountain.meadows", "Reed", 
                                     "Swamps","Desert","Freshwater",
                                     "Marine")],1,max)

### avonnet data ####s

avonet <- read.csv("~/Dropbox/DOF/traits/Avonet/Supplementary_dataset_1.csv", 
                   as.is = TRUE, sep=";")
species[!species %in% avonet$Species1]

avonet$Species1[avonet$Species1 %in% c("Acanthis flammea")] <- "Acanthis flammea/A. cabaret"
avonet$Species1[avonet$Species1 %in% c("Larus ridibundus")] <- "Chroicocephalus ridibundus"
avonet$Species1[avonet$Species1 %in% c("Corvus monedula")] <- "Coloeus monedula"
avonet$Species1[avonet$Species1 %in% c("Sylvia communis")] <- "Curruca communis"
avonet$Species1[avonet$Species1 %in% c("Sylvia curruca")] <- "Curruca curruca"
avonet$Species1[avonet$Species1 %in% c("Saxicola torquatus")] <- "Saxicola rubicola"

#add on old sub species for crows
avonet2 <-  subset(avonet, Species1=="Corvus corone")
avonet2$Species1 <- "Corvus cornix"
avonet <- rbind(avonet, avonet2)

#clean 
avonet <- avonet %>%
            rename(Species = "Species1") %>%
            filter(Species %in% species) %>%
            select(Species, Family1, Order1, 
                   Sequence, Avibase.ID1,
                   Mass, Habitat, Habitat.Density,
                   Migration, Trophic.Level, Trophic.Niche,
                   Primary.Lifestyle, Range.Size, Centroid.Latitude)

write.table(avonet,file="avonet_traits.txt",sep="\t",row.names=FALSE)

alltraits <- inner_join(alltraits, avonet, by="Species")

### save #####  

saveRDS(alltraits, "traits.rds")

### end ####
