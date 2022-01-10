library(raster)
library(terra)
library(tidyverse)
library(sf)
library(unmarked)
library(formula.tools)

Puntos_Hull <- read_csv("https://raw.github.com/derek-corcoran-barrios/LayerCreationBuffer/main/Coords.csv") %>% 
  mutate(geometry = str_remove_all(str_remove_all(str_remove_all(geometry, "c"), "\\("), "\\)")) 

Puntos_Hull$Lon <-  str_split(Puntos_Hull$geometry, pattern = ",", simplify = T)[,1] %>% as.numeric()
Puntos_Hull$Lat <-  str_split(Puntos_Hull$geometry, pattern = ", ", simplify = T)[,2] %>% as.numeric()

Puntos_Hull <- Puntos_Hull %>% 
  dplyr::select(-geometry) %>% 
  st_as_sf(coords = c(2,3), crs = 4326) %>% 
  st_transform(crs = "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs") %>%
  terra::vect()

Distancias <- round(seq(from = 30, to = 5000, length.out = 10), -2)
Distancias[1] <- 30

DistanciaRio <- readRDS("DistanceToRiver.RDS") %>% 
  raster::projectRaster(crs = "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs") %>% 
  terra::rast()
Altura <- readRDS("Alt.rds") %>% 
  raster::projectRaster(crs = "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs") %>% 
  terra::rast()


Layers <- list()
OccuVars <- list()

for(i in 1:length(Distancias)){
  Layers[[i]] <- terra::rast(paste0("/vsicurl/https://raw.github.com/derek-corcoran-barrios/LayerCreationBuffer/main/Proportions_", Distancias[i],".tif")) %>% 
    terra::resample(Altura)
  Layers[[i]] <- c(Layers[[i]], Altura, DistanciaRio)
  names(Layers[[i]])[c(10,11)] <- c("altura", "distancia_rio")
  OccuVars[[i]] <- terra::extract(Layers[[i]], Puntos_Hull)
  names(Layers[[i]])[c(10,11)]
}


data_reg_Inv <- read_rds("Occdata_regInv.rds") %>% 
  dplyr::select(matches("Pelecanus_thagus|Larus_dominicanus|Coragyps_atratus|Larosterna_inca|Columba_livia|Sephanoides_sephaniodes"))


Species_Names <- unique(gsub('[0-9]+', '', colnames(data_reg_Inv)))

N_Species <- length(Species_Names)

data_det_Inv <-read_rds("Occdata_detInv.rds")

batchoccu2 <- function(pres, sitecov, obscov, spp, form, SppNames = NULL, dredge = FALSE) {
  if(is.null(SppNames)){
    SppNames <- paste("species", 1:spp, sep =".")
  }
  secuencia <- c(1:spp)*(ncol(pres)/spp)
  secuencia2<-secuencia-(secuencia[1]-1)
  models <- vector('list', spp)
  fit <- matrix(NA, nrow(pres), spp)
  Mods <- list()
  if(is.null(SppNames)){
    colnames(fit) <- paste("species", 1:spp, sep =".")
  }else if(class(SppNames) == "character"){
    colnames(fit) <- SppNames
  }
  if (dredge == FALSE) {
    for(i in 1:length(secuencia)) {
      data <- pres[, secuencia2[i]:secuencia[i]]
      data2 <- unmarkedFrameOccu(y = data, siteCovs = sitecov, obsCovs = obscov)
      try({
        models[[i]] <- occu(as.formula(form), data2)
      }, silent = T)
      try({
        fit[, i] <- suppressWarnings(predict(models[[i]], type = "state", newdata = sitecov))$Predicted
      }, silent = T)
      Mods = NULL
      print(paste("Species", as.character(i), "ready!"))
    }
  }
  else {
    for(i in 1:length(secuencia)) {
      data <- pres[, secuencia2[i]:secuencia[i]]
      data2 <- unmarkedFrameOccu(y = data, siteCovs = sitecov, obsCovs = obscov)
      try({
        #Partimos en dos Detección y occupancia
        form <- as.character(form)
        Div <- str_squish(form) %>% str_remove_all(" ")  %>% stringr::str_split(pattern = "~", simplify = T)
        
        ### Separamos dos formulas Occupancia y Deteccion
        
        Det <- Div[length(Div) - 1]
        
        VarDet <- str_split(Det, "\\+", simplify = T) %>% as.character()
        
        Fs <- list()
        
        print(paste("Starting to fit detection models for species", i, "of", length(secuencia)))
        
        for(x in 1:(length(VarDet) + 1)){
          if(x == (length(VarDet) + 1)){
            Formulas <- data.frame(Form = "~1 ~ 1", AICc = NA)
            Formulas$AICc[j] <- try(MuMIn::AICc(occu(as.formula("~1 ~1"), data2)), silent = T)
          }else{
            Test <- combn(VarDet, x, simplify = F)
            Formulas <- data.frame(Form = rep(NA, length(Test)), AICc = rep(NA, length(Test)))
            for(j in 1:length(Test)){
              Temp <- paste("~", paste(Test[[j]], collapse = " + "), "~ 1") 
              Formulas$Form[j] <- Temp
              Temp <- as.formula(Temp)
              Formulas$AICc[j] <- try(MuMIn::AICc(occu(Temp, data2)), silent = T) 
              gc()
            }
          }
          
          Fs[[x]] <- suppressWarnings(Formulas %>% mutate(AICc = as.numeric(AICc)) %>% dplyr::filter(!is.na(AICc)) %>% arrange(AICc))
          message(paste("finished for", x, "number of variables"))
        }
        
        Fs <- suppressWarnings(purrr::reduce(Fs, bind_rows) %>% arrange(AICc))
        
        Selected <- Fs$Form[1] %>% str_split("~", simplify = T) %>% as.character()
        Selected <- Selected[length(Selected) - 1] %>% str_squish()
        
        print(paste("Detection model for species", i, "is", Selected))
        
        Occup <- Div[length(Div)]
        
        VarOccup <- str_split(Occup, "\\+", simplify = T) %>% as.character()
        
        Fs <- list()
        
        print(paste("Starting to fit occupancy models for species", i, "of", length(secuencia)))
        
        for(x in 1:(length(VarOccup) + 1)){
          if(x == (length(VarOccup) + 1)){
            Formulas <- data.frame(Form = paste("~",Selected, "~ 1"), AICc = NA)
            Formulas$AICc[j] <- try(MuMIn::AICc(occu(as.formula(paste("~",Selected, "~ 1")), data2)), silent = T)
          }else{
            Test <- combn(VarOccup, x, simplify = F)
            Formulas <- data.frame(Form = rep(NA, length(Test)), AICc = rep(NA, length(Test)))
            for(j in 1:length(Test)){
              Temp <- paste("~", Selected, "~", paste(Test[[j]], collapse = " + ")) 
              Formulas$Form[j] <- Temp
              Temp <- as.formula(Temp)
              Formulas$AICc[j] <- try(MuMIn::AICc(occu(Temp, data2)), silent = T) 
              if((j %% 100) == 0){
                message(paste(j, "of", length(Test), "Ready"))
                gc()
              }
            }
          }
          
          Fs[[x]] <- suppressWarnings(Formulas %>% mutate(AICc = as.numeric(AICc)) %>% dplyr::filter(!is.na(AICc)) %>% arrange(AICc))
          message(paste("finished for", x, "number of variables", Sys.time()))
        }
        
        Fs <- suppressWarnings(purrr::reduce(Fs, bind_rows) %>% arrange(AICc))
        
        Mods[[i]] <- Fs
        
        
        Best <- Fs$Form[1]
        
        models[[i]] <- occu(as.formula(Best), data2)
        #dredged <- suppressWarnings(dredge(occu(form, data2)))
        # select the first model and evaluate
        #models[[i]] <- eval(getCall(dredged, 1))
        
      }, silent = T)
      try({
        #predictions for the best model
        fit[, i] <- suppressWarnings(predict(models[[i]], type = "state", newdata = sitecov))$Predicted
      }, silent = T)
      
      print(paste("Species", as.character(i), "ready!"))
    }
  }
  if(is.null(SppNames)){
    names(models) <- paste("species", 1:spp, sep =".")
  }else if(class(SppNames) == "character"){
    names(models) <- SppNames
  }
  
  if(is.null(SppNames)){
    names(Mods) <- paste("species", 1:spp, sep =".")
  }else if(class(SppNames) == "character" & !is.null(Mods)){
    names(Mods) <- SppNames
  }
  
  cond <- sapply(models, function(x) !is.null(x))
  models <- models[cond]
  fit <- fit[,cond]
  Not <- SppNames[!(cond)]
  if(sum(!cond) >= 1){
    message(paste("species", paste(Not, collapse = ", "), "did not converge, try with less variables"))
  }
  result <- list(Covs = sitecov, models = models, fit = fit, Mods = Mods)
  class(result)<- "batchoccupancy"
  return(result)
}


######## Con solo algunas variables


Results_Inv <- list()

for(i in 1:length(Distancias)){
  print(paste("Starting distance", Distancias[i], Sys.time()))
  Results_Inv[[i]] <- batchoccu2(pres = data_reg_Inv, sitecov = OccuVars[[i]], obscov = data_det_Inv, 
                             spp=N_Species, 
                             form= "~ Temperatura + Humedad + DirViento + RapViento + Agua ~ bosque_nativo + cultivos + grava + oceano + pastizales + matorrales + sup_impermeables + suelo_arenoso + plantacion_de_arboles + altura + distancia_rio", 
                             dredge=TRUE,  
                             SppNames = Species_Names)
}


Table_Inv <- list()
Template <- data.frame(Species = NA, Distance = NA, Model = NA, AICc = NA)

for(i in 1:length(Results_Inv)){
  TempDistance <- list()
  for(j in 1:length(Results_Inv[[i]]$models)){
    Temp <- Template
    Temp$Species <- (names(Results_Inv[[i]]$models))[j] 
    Temp$Distance <- Distancias[i]
    Temp$Model <- Results_Inv[[i]]$models[[j]]@formula %>% as.character()
    Temp$AICc <- Results_Inv[[i]]$models[[j]] %>% MuMIn::AICc()
    TempDistance[[j]] <- Temp
  }
  Table_Inv[[i]] <- TempDistance %>% purrr::reduce(bind_rows)
}

Table_Inv <- Table_Inv %>% 
  purrr::reduce(bind_rows) %>% 
  group_split(Species) %>% 
  purrr::map(~mutate(.x, delta_AICc = AICc - min(AICc))) %>% 
  purrr::map(~arrange(.x, delta_AICc)) %>% 
  purrr::map(~dplyr::filter(.x, delta_AICc <= 2)) %>% 
  purrr::map(~dplyr::filter(.x, Distance == min(Distance))) %>% 
  purrr::reduce(bind_rows) 

## PRedicciones invierno

dir.create("Results_Winter")
dir.create("Results_Winter/Tables")
dir.create("Results_Winter/Plots")
dir.create("Results_Winter/Rasters")

write_csv(Table_Inv, "Results_Winter/Tables/Table_Models_inv.csv")

Stack_Species_Inv <- list()

for(i in 1:nrow(Table_Inv)){
  Cond <- (Distancias == Table_Inv$Distance[i])
  Index <- (1:length(Distancias))[Cond]
  Stack_Species_Inv[[i]] <- predict(Results_Inv[[Index]]$models[[i]] ,raster::stack(Layers[[Index]]), type = "state")[[1]]
  raster::writeRaster(Stack_Species_Inv[[i]], paste0("Results_Winter/Rasters/",Table_Inv$Species[i],".tif"))
  message(paste(i, "of",nrow(Table_Inv), Sys.time()))
}

MaskRaster <- raster(Altura)
MaskRaster[!is.na(MaskRaster)] <- 1

for(i in 1:length(Stack_Species_Inv)){
  Stack_Species_Inv[[i]] <- Stack_Species_Inv[[i]]*MaskRaster
}

Stack_Species_Inv <- Stack_Species_Inv %>% purrr::reduce(stack)

names(Stack_Species_Inv) <- Table_Inv$Species

png("Results_Winter/Plots/Species.png",res = 300, width = 2000, height = 2000)
plot(Stack_Species_Inv, colNA = "black")
dev.off()
### Primavera

data_reg_Prim <- read_rds("Occdata_regPRIM.rds") %>% 
  dplyr::select(matches("Leucophaeus_pipixcan|Larus_dominicanus|Columba_livia|Larosterna_inca|Phalacrocorax_bougainvillii"))

Species_Names <- unique(gsub('[0-9]+', '', colnames(data_reg_Prim)))

N_Species <- length(Species_Names)

data_det_Prim <-read_rds("Occdata_detPrim.rds")

Results_Prim <- list()

for(i in 1:length(Distancias)){
  print(paste("Starting distance", Distancias[i], Sys.time()))
  Results_Prim[[i]] <- batchoccu2(pres = data_reg_Prim, sitecov = OccuVars[[i]], obscov = data_det_Prim, 
                                 spp=N_Species, 
                                 form= "~ Temperatura + Humedad + DirViento + RapViento + Agua ~ bosque_nativo + cultivos + grava + oceano + pastizales + matorrales + sup_impermeables + suelo_arenoso + plantacion_de_arboles + altura + distancia_rio", 
                                 dredge=TRUE,  
                                 SppNames = Species_Names)
}


Table_Prim <- list()
Template <- data.frame(Species = NA, Distance = NA, Model = NA, AICc = NA)

for(i in 1:length(Results_Prim)){
  TempDistance <- list()
  for(j in 1:length(Results_Prim[[i]]$models)){
    Temp <- Template
    Temp$Species <- (names(Results_Prim[[i]]$models))[j] 
    Temp$Distance <- Distancias[i]
    Temp$Model <- Results_Prim[[i]]$models[[j]]@formula %>% as.character()
    Temp$AICc <- Results_Prim[[i]]$models[[j]] %>% MuMIn::AICc()
    TempDistance[[j]] <- Temp
  }
  Table_Prim[[i]] <- TempDistance %>% purrr::reduce(bind_rows)
}

Table_Prim <- Table_Prim %>% 
  purrr::reduce(bind_rows) %>% 
  group_split(Species) %>% 
  purrr::map(~mutate(.x, delta_AICc = AICc - min(AICc))) %>% 
  purrr::map(~arrange(.x, delta_AICc)) %>% 
  purrr::map(~dplyr::filter(.x, delta_AICc <= 2)) %>% 
  purrr::map(~dplyr::filter(.x, Distance == min(Distance))) %>% 
  purrr::reduce(bind_rows)


### Final models and projections
## PRedicciones Primavera

dir.create("Results_Spring")
dir.create("Results_Spring/Tables")
dir.create("Results_Spring/Plots")
dir.create("Results_Spring/Rasters")

write_csv(Table_Prim, "Results_Spring/Tables/Table_Models_prim.csv")

Stack_Species_Prim <- list()

for(i in 1:nrow(Table_Prim)){
  Cond <- (Distancias == Table_Prim$Distance[i])
  Index <- (1:length(Distancias))[Cond]
  Stack_Species_Prim[[i]] <- predict(Results_Prim[[Index]]$models[[i]] ,raster::stack(Layers[[Index]]), type = "state")[[1]]
  raster::writeRaster(Stack_Species_Prim[[i]], paste0("Results_Spring/Rasters/",Table_Prim$Species[i],".tif"))
  message(paste(i, "of",nrow(Table_Prim), Sys.time()))
}

for(i in 1:length(Stack_Species_Prim)){
  Stack_Species_Prim[[i]] <- Stack_Species_Prim[[i]]*MaskRaster
}

Stack_Species_Prim <- Stack_Species_Prim %>% purrr::reduce(stack)

names(Stack_Species_Prim) <- Table_Prim$Species

png("Results_Spring/Plots/Species.png",res = 300, width = 2000, height = 2000)
plot(Stack_Species_Prim, colNA = "black")
dev.off()