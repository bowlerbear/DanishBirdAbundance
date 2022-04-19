# script to get a phylogenetic consensus tree to use for phylogenetic analyses

library(tidyverse)
library(readr)
library(phylolm)
library(phylosignal)
library(phylobase)
library(phytools)

### get species from danish study ####
species <- read_csv("DOF_species.txt")


# read in birdtree taxonomy of avonet
# this is from supplementary dataset 1
# and I saved the AVONET3_BirdTree tab out separately
# because this is what we want to match to
avonet_birdtree <- read_csv("avonet_birdtree.csv")

# now need to match the taxonomy
matched <- species %>%
  left_join(., avonet_birdtree %>%
              rename(Species=Species3)) %>%
  dplyr::filter(complete.cases(Family3))

# 103 match exactly
# 20 do not match
not_matched <- species %>%
  dplyr::filter(!Species %in% matched$Species)

taxonomy_fix <- not_matched %>%
  mutate(species_birdtree=c("Carduelis flammea", "Carduelis chloris", 
                            "Larus ridibundus", "Corvus monedula", 
                            "Corvus corone", "Sylvia communis",
                            "Sylvia curruca", "Parus caeruleus",
                            "Miliaria calandra", "Larus canus",
                            "Carduelis cannabina", "Parus cristatus",
                            "Anas strepera", "Parus ater", "Parus montanus",
                            "Parus palustris", "Saxicola torquatus",
                            "Anas clypeata", "Anas querquedula",
                            "Carduelis spinus"))

not_matched_fixed <- taxonomy_fix %>%
  rename(Species3=species_birdtree) %>%
  left_join(., avonet_birdtree) 

# now combine the two together into 1 dataset
# a bit hacky for sure
all_combined <- matched %>%
  bind_rows(not_matched_fixed) %>%
  rename(danish_species=Species) %>%
  mutate(Species3=ifelse(is.na(Species3)==TRUE, danish_species, Species3)) %>%
  dplyr::select(1, 37, 2:36) %>%
  rename(TipLabel=Species3)

# does TipLabel match the number of rows?
length(unique(all_combined$TipLabel))==nrow(all_combined)  

# NO!
# which ones are repeating?
all_combined %>%
  group_by(TipLabel) %>%
  summarize(N=n()) %>%
  arrange(desc(N))

# now will just slice one row for each of these two TipLabels that
# have mutliple measures/species corresponding in danish dataset
all_combined <- all_combined %>%
  group_by(TipLabel) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(TipLabel=gsub(" ", "_", TipLabel))
  
# does TipLabel match the number of rows?
length(unique(all_combined$TipLabel))==nrow(all_combined)  

# let's write the file out
saveRDS(all_combined, "avonet_traits_birdtree_taxonomy.RDS")

# YES! Let's get a phylogenetic tree!
#########################################
##########################################
#### Now start phylo stuff ###############

# function to read one tree in
read_one_tree <- function(path, x=1){
  
  one_bird_tree <- ape::read.tree(file = "phylo/phy.tre")[[x]]
  
  return(one_bird_tree)
}

bird_tree <- read_one_tree()


# function to read all trees in
read_all_trees <- function(path){
  
  ape::read.tree(file = "phylo/phy.tre")
  
}

all_trees <- read_all_trees()

# a function to subset the tree to the tips of the species described above
subset_tree <- function(bird_tree, dataset) {
  
  non_danish_sp <- bird_tree$tip.label[!bird_tree$tip.label %in% dataset$TipLabel]
  
  danish_bird_tree <- drop.tip(bird_tree, non_danish_sp)
  
  return(danish_bird_tree)
}


# now we create a phylogenetic tree for the species that we have
danish_tree <- subset_tree(bird_tree, all_combined)
plot(danish_tree)

# need to get a consensus tree to run the phylogenetic analyses on
# first subset all trees to the 603 species
non_danish_sp <- bird_tree$tip.label[!bird_tree$tip.label %in% all_combined$TipLabel]

subset_trees <- lapply(all_trees, drop.tip, tip=non_danish_sp)

# now make a consensus tree
# this takes a while to run...
con_tree <- consensus.edges(subset_trees, consensus.tree=consensus(subset_trees, p=0.5, check.labels=TRUE))

saveRDS(con_tree, "consensus_tree.RDS")

con_tree <- readRDS("consensus_tree.RDS")

##########################
##########################
####### example model

# will just use the "Total.individuals" as a fake response variable
ex_mod_dat <- all_combined %>%
  dplyr::select(TipLabel, Total.individuals, Mass)

# make tip label rownnames
row.names(ex_mod_dat) <- ex_mod_dat$TipLabel

# fit a simple linear model 
mod1 <- lm(Total.individuals ~ log10(Mass), data=ex_mod_dat)
summary(mod1)

# phylogenetic model version
phy_mod1 <- phylolm(Total.individuals ~ log10(Mass), data=ex_mod_dat,
                    phy=con_tree, na.action="na.fail")
summary(phy_mod1)

# switches from being super strong effect in non phylo model to no effect in phylo model.
# as expected since mass has strong phylogenetic dependence




##################################################
##################################################
##################################################
############### Get color data
##################################################
















