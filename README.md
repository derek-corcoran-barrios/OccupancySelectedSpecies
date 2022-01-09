Occupancy models for selected species
================
Giorgia Graells and Derek Corcoran
2022-01-09

# 1 Objective

This repository serves to document and store the analyses and results
for the manuscript **Exploring habitat use of terrestrial and marine
birds in urban coastal areas** sent to the jounrnal **Frontiers Ecology
And Evolution**

## 1.1 Methods

First we load the required packages

``` r
# For manipulating and reading raster datasets
library(raster)
library(terra)
# For cleaning datasets
library(tidyverse)
# For managing vector spatial datasets
library(sf)
# For caclulating occupancy models
library(unmarked)
# For selecting models
library(MuMIn)
```

Then we load the coordinates of the sampling sites for the surveys and
transform them in to a `SpatVector` object:

``` r
Puntos_Hull <- read_csv("https://raw.github.com/derek-corcoran-barrios/LayerCreationBuffer/main/Coords.csv") %>% 
  mutate(geometry = str_remove_all(str_remove_all(str_remove_all(geometry, "c"), "\\("), "\\)")) 

Puntos_Hull$Lon <-  str_split(Puntos_Hull$geometry, pattern = ",", simplify = T)[,1] %>% as.numeric()
Puntos_Hull$Lat <-  str_split(Puntos_Hull$geometry, pattern = ", ", simplify = T)[,2] %>% as.numeric()

Puntos_Hull <- Puntos_Hull %>% 
  dplyr::select(-geometry) %>% 
  st_as_sf(coords = c(2,3), crs = 4326) %>% 
  st_transform(crs = "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs") %>%
  terra::vect()
```

Then we generate a vector of the distances used to calculate the
proportion of landuse in meters as seen in the LayerCreationBuffer
repository
