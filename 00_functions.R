addCoords <- function(df){
  
  require(sf)
  squares32 <- st_read(dsn = "data/TTT", layer = "transect squares utm32")
  squares33 <- st_read(dsn = "data/TTT", layer = "transect squares utm33") %>%
    st_transform(.,st_crs(squares32))
  squares <- bind_rows(squares32, squares33) %>%
    filter(kvadratnr %in% data$kvadratnr) %>%
    filter(!duplicated(kvadratnr))
  
  coordsDF <- st_coordinates(st_centroid(squares)) %>%
              as_tibble() %>%
              add_column(kvadratnr = squares$kvadratnr)
 
  df <- inner_join(df,coordsDF, by="kvadratnr")
   return(df)
}

getDistanceData <- function(myspecies){
  
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
  
  #get number of birds seen at each distance
  allData %>% 
    dplyr::select(Year, X.0,X.1,X.2) %>%
    pivot_longer(starts_with("X."),names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = "25", X.1 = "50", X.2 = "100")) %>%
    dplyr::group_by(Year, distance) %>%
    dplyr::summarise(total = sum(nu)) %>%
    dplyr::mutate(distance = factor(distance, levels = c("25", "50", "100"))) %>%
    ungroup() %>%
    add_column(Species = myspecies)
  
}

fitModel <- function(myspecies){
  
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
  
  ### format data for unmarked ####
  
  allData$lines_pathroad2 <- allData$lines_pathroad^2
  covariates <- data.frame(scale(allData[,c("skydaekke","regn","vind","Year",
                                            "lines_path","lines_road","lines_forest",
                                            "lines_pathroad","lines_pathroad2",
                                            "squares_forest","squares_agri_int",
                                            "squares_urban", "squares_agri_ext",
                                            "squares_freshwater")]))
  #covariates$fYear <- factor(allData$Year)
  
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
  
  #null model
  fm0 <- distsamp(~ 1
                  ~1, 
                  data = unmarkDF, 
                  keyfun = "halfnorm",
                  output = "density",
                  unitsOut = "kmsq")
  
  #with covariates
  fm1 <- distsamp(~ lines_path + lines_road
                  ~ squares_forest + squares_agri_int + squares_urban + 
                    squares_agri_ext + squares_freshwater, 
                  data = unmarkDF, 
                  keyfun = "halfnorm",
                  output = "density",
                  unitsOut = "kmsq")
  
  #null state model
  fm2 <- distsamp(~ lines_path + lines_road 
                  ~ 1, 
                  data = unmarkDF, 
                  keyfun = "halfnorm",
                  output = "density",
                  unitsOut = "kmsq")
  
  
  #null detection model
  fm3 <- distsamp(~ 1
                  ~ squares_forest + squares_agri_int + squares_urban + 
                    squares_agri_ext + squares_freshwater, 
                  data = unmarkDF, 
                  keyfun = "halfnorm",
                  output = "density",
                  unitsOut = "kmsq")
  
  list(species = myspecies,
        fm0 = fm0,
        fm1 = fm1, 
        fm2 = fm2,
        fm3 = fm3)
  
}

extractModel <- function(myModels,modeltype, output="ESW", sep=TRUE){
  
  
  #and just one species
  dataS <- data %>% 
    filter(Species == myModels[["species"]]) %>%
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
  
  ### format data for unmarked ####
  
  allData$lines_pathroad2 <- allData$lines_pathroad^2
  covariates <- data.frame(scale(allData[,c("skydaekke","regn","vind","Year",
                                            "lines_path","lines_road","lines_forest",
                                            "lines_pathroad","lines_pathroad2",
                                            "squares_forest","squares_agri_int",
                                            "squares_urban", "squares_agri_ext",
                                            "squares_freshwater")]))
  #covariates$fYear <- as.factor(allData$Year)
  
  #distance data
  temp <- allData[,c("X.0","X.1","X.2")]
  temp[is.na(temp)] <- 0
  colSums(temp)
  
  #select model
  
  model <- myModels[[modeltype]]
  
  #decide on output
  if(output=="state"){
    
    stateDF <- data.frame(Species = myModels[["species"]],
                          param = names(coef(model,type='state')),
                          coef = as.numeric(coef(model,type='state')),
                          coef_se = as.numeric(SE(model,type='state')))
    return(stateDF)
    
  }else if (output=="detection"){
    
    detectionDF <- data.frame(Species = myModels[["species"]],
                              param = names(coef(model,type='det')),
                              coef = as.numeric(coef(model,type='det')),
                              coef_se = as.numeric(SE(model,type='det')))
    
    return(detectionDF)
    
    #use model to predict sigma
  }else if(output=="ESW"){
    
    #newData = data.frame(lines_pathroad = min(covariates$lines_pathroad))
    if(sep==TRUE){
    newData = data.frame(lines_path = min(covariates$lines_path),
                          lines_road = min(covariates$lines_road))
    }else{
      newData = data.frame(lines_pathroad = min(covariates$lines_pathroad))
    }
    
    log_sigma <- predict(model, type="det", newdata = newData)
    
    #get effective strip width
    log_sigma$ESW <- sapply(log_sigma$Predicted, function(x){
      integrate(gxhn, 0, 101, sigma = x)$value
    }) 
    #head(log_sigma)#sigma is 38.2, ESW is 47 which makes sense
    #backTransform(fm1, type="det")  #same as what is in predict             
    #sqrt((pi * 38.51699^2)/2)
    
    #add on species
    log_sigma$Species <- myModels[["species"]]
    
    return(log_sigma)
    
  }else if(output=="density"){
    
    densityDF <- predict(model, type="state", newdata = covariates)
    densityDF <- bind_cols(densityDF, covariates)
    densityDF$kvadratnr <- allData$kvadratnr
    
    #add on species
    densityDF$Species <- myModels[["species"]]
    return(densityDF)
    
  } 
  
}

