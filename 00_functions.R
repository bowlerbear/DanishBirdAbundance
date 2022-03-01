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
  allData <- left_join(info, dataS)
  
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
  allData <- left_join(info, dataS)
  
  #add on environmental data
  allData <- inner_join(allData, environData)
  
  ### format data for unmarked ####
  
  covariates <- data.frame(scale(allData[,c("skydaekke","regn","vind","Year",
                                            "lines_path","lines_road","lines_forest",
                                            "lines_pathroad",
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
  fm1 <- distsamp(~ fYear + lines_pathroad
                  ~ fYear + squares_forest + squares_agri_int + squares_urban + 
                    squares_agri_ext + squares_freshwater, 
                  data = unmarkDF, 
                  keyfun = "halfnorm",
                  output = "density",
                  unitsOut = "kmsq")
  
  #null model
  fm0 <- distsamp(~ 1
                  ~1, 
                  data = unmarkDF, 
                  keyfun = "halfnorm",
                  output = "density",
                  unitsOut = "kmsq")
  
  
  #null state model
  fm2 <- distsamp(~ fYear + lines_pathroad
                  ~1, 
                  data = unmarkDF, 
                  keyfun = "halfnorm",
                  output = "density",
                  unitsOut = "kmsq")
  
  list(fm0 = fm0,
       fm1 = fm1, 
       fm2 = fm2)
  
}

extractModel <- function(myModels,modeltype, output="ESW"){
  
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
  allData <- left_join(info, dataS)
  
  #add on environmental data
  allData <- inner_join(allData, environData)
  
  ### format data ####
  
  covariates <- data.frame(scale(allData[,c("skydaekke","regn","vind","Year",
                                            "lines_path","lines_road","lines_forest",
                                            "lines_pathroad",
                                            "squares_forest","squares_agri_int",
                                            "squares_urban", "squares_agri_ext",
                                            "squares_freshwater")]))
  covariates$fYear <- as.factor(allData$Year)
  
  #distance data
  temp <- allData[,c("X.0","X.1","X.2")]
  temp[is.na(temp)] <- 0
  colSums(temp)
  
  #select model
  
  model <- myModels[[modeltype]]
  
  #decide on output
  if(output=="state"){
    
    stateDF <- data.frame(Species = myspecies,
                          param = names(coef(model,type='state')),
                          coef = as.numeric(coef(model,type='state')),
                          coef_se = as.numeric(SE(model,type='state')))
    return(stateDF)
    
  }else if (output=="detection"){
    
    detectionDF <- data.frame(Species = myspecies,
                              param = names(coef(model,type='det')),
                              coef = as.numeric(coef(model,type='det')),
                              coef_se = as.numeric(SE(model,type='det')))
    
    return(detectionDF)
    
    #use model to predict sigma
  }else if(output=="ESW"){
    
    newData = expand_grid(fYear = sort(unique(allData$Year)),
                          lines_pathroad = min(covariates$lines_pathroad)) %>%
      mutate(fYear = as.factor(fYear)) %>%
      as.data.frame()
    
    log_sigma <- predict(model, type="det", newdata = newData)
    
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
    
  }else if(output=="density"){
    
    densityDF <- predict(model, type="state", newdata = covariates)
    
    #add on species
    densityDF$Species <- myspecies
    return(densityDF)
    
  } 
  
}

