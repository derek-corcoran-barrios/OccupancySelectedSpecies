library(raster)
library(terra)
library(tidyverse)
library(sf)

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

Layers <- list()
OccuVars <- list()

for(i in 1:length(Distancias)){
  Layers[[i]] <- terra::rast(paste0("/vsicurl/https://raw.github.com/derek-corcoran-barrios/LayerCreationBuffer/main/Proportions_", Distancias[i],".tif"))
  OccuVars[[i]] <- terra::extract(Layers[[i]], Puntos_Hull)
}


data_reg_Inv <- read_rds("Occdata_regInv.rds") %>% 
  dplyr::select(matches("Pelecanus_thagus|Larus_dominicanus|Coragyps_atratus|Larosterna_inca|Columba_livia|Sephanoides_sephaniodes"))


data_reg_Prim <- read_rds("Occdata_regPRIM.rds") %>% 
  dplyr::select(matches("Leucophaeus_pipixcan|Larus_dominicanus|Columba_livia|Larosterna_inca|Phalacrocorax_bougainvillii"))

Species_Names <- unique(gsub('[0-9]+', '', colnames(data_reg_Inv)))

N_Species <- length(Species_Names)

data_det <-read_rds("Occdata_detInv.rds")

data_ocu <-read_rds("Occdata_occu.rds")
data_ocu <- data_ocu %>% select(-Sitio)
colnames(data_ocu) <- str_replace_all(colnames(data_ocu), pattern = " ", "_")

data_ocu2 <- data_ocu %>% mutate_if(is.numeric, scale) #centro y escalamiento de variables-media 0, varianza1

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

# BUFFER 30 M DE LOS SITIOS DE MUESTREO
OccuInv37_30 <- batchoccu2(pres = data_reg, sitecov = data_ocu, obscov = OccuVars[[i]], 
                           spp=N_Species, 
                           form= "~ Temperatura + Humedad + DirViento + RapViento + Agua ~ pastizales + matorrales+ sup_impermeables + oceano + suelo_arenoso + Buffer_30_Rocas + grava", 
                           dredge=TRUE,  
                           SppNames = Species_Names)

OccuInv37_600 <- batchoccu2(pres = data_reg, 
                            sitecov = data_ocu, 
                            obscov = data_det, 
                            spp=N_Species , 
                            form= "~ Temperatura + Humedad + DirViento + RapViento + Agua ~ bosque_nativo + pastizales + matorrales+ Buffer_600_Humedales + Buffer_600_Sup_impermeables+ Buffer_600_Oceano+ Buffer_600_Suelo_arenoso+Buffer_600_Rocas+ Buffer_600_Grava", 
                            dredge=TRUE,  
                            SppNames = Species_Names)