## Functions used to in the .inputObjects event:
# prepNTEMSDominantSpecies
# prepSoilTexture
# prepNdeposition
# rvestAlbedoTable
# prepCo2Concentration
# prepEPC
# prepClimate

# Extract dominant species from NTEMS layers
prepNTEMSDominantSpecies <- function(year, destinationPath, cropTo, projectTo, maskTo, method = "mode"){
  # prepare url and targetFile for prepInput call
  domSppURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_Tree_Species_Classification_", year, ".zip")
  domSppTF <- paste0("CA_Forest_Tree_Species_", year, ".tif")  
  
  # prepInput
  dominantSpecies <- prepInputs(
    targetFile = domSppTF,
    url = domSppURL,
    destinationPath = destinationPath,
    cropTo = cropTo,
    projectTo = projectTo,
    maskTo = maskTo,
    method = method,
    overwrite = TRUE
  )
  # set 0 to NA
  NAflag(dominantSpecies) <- 0
  
  # Transform NTEMS code to speciesCode
  cls <- unique(LandR::sppEquivalencies_CA[, c("NTEMS_Species_Code", "LandR")]) |> na.omit()
  names(cls) <- c("id", "category")
  levels(dominantSpecies) <- cls
  
  return(dominantSpecies)
}

# Extract % of sand, % of clay and % of silt from CanSIS dataset
# Takes the average across layers 0-30cm
prepSoilTexture <- function(destinationPath, to){
  ## Sand %
  # layer 0-5 cm
  sand0_5 <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Sand/Sand_X0_5_cm_100m1980-2000v1.tif",
    targetFile = "Sand_X0_5_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to,
    overwrite = TRUE
  ) |> Cache()
  sand5_15 <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Sand/Sand_X5_15_cm_100m1980-2000v1.tif",
    targetFile = "Sand_X5_15_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to,
    overwrite = TRUE
  ) |> Cache()
  sand15_30 <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Sand/Sand_X15_30_cm_100m1980-2000v1.tif",
    targetFile = "Sand_X15_30_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to,
    overwrite = TRUE
  ) |> Cache()
  
  sand <- round((5/30) * sand0_5 + (10/30) * sand5_15 + (15/30) * sand15_30, digit = -1)
  
  clay0_5 <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Clay/Clay_X0_5_cm_100m1980-2000v1.tif",
    targetFile = "Clay_X0_5_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to,
    overwrite = TRUE
  ) |> Cache()
  
  clay5_15 <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Clay/Clay_X5_15_cm_100m1980-2000v1.tif",
    targetFile = "Clay_X5_15_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to,
    overwrite = TRUE
  ) |> Cache()
  
  clay15_30 <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Clay/Clay_X15_30_cm_100m1980-2000v1.tif",
    targetFile = "Clay_X15_30_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to,
    overwrite = TRUE
  ) |> Cache()
  
  clay <- round((5/30) * clay0_5 + (10/30) * clay5_15 + (15/30) * clay15_30, digit = -1)
  
  silt <- 100 - (sand + clay)
  soilTexture <- c(sand, silt, clay)
  names(soilTexture) <- c("sand", "silt", "clay")
  return(soilTexture)
}

# Extract N-deposition data for 2 years.
prepNdeposition <- function(destinationPath, to, year1, year2){
  Ndeposition1 <- prepInputs(
    targetFile = paste0("mean_totN_", year1, "_hm.tif"),
    archive = "Global_N_deposition_grid_dataset_2008_2020.rar",
    overwrite = TRUE,
    url = "https://springernature.figshare.com/ndownloader/files/48644623",
    cropTo = to,
    projectTo = to,
    destinationPath = destinationPath,
    fun = "terra::rast"
  )
  # reprojecting can create holes
  w <- sum(is.na(values(Ndeposition1)))
  while (w != 0){
    w <- round(sqrt(w))
    # make sure w is a odd number
    if(w %% 2 != 1){
      w = w+1
    }
    Ndeposition1 <- focal(Ndeposition1, w = w, fun = 'mean', na.policy = 'only')
    w <- sum(is.na(values(Ndeposition1)))
  }
  # mask to study area
  Ndeposition1 <- maskTo(Ndeposition1, to) |> round()
  
  # N deposition for the second time step
  Ndeposition2 <- prepInputs(
    targetFile = paste0("mean_totN_", year2, "_hm.tif"),
    archive = "Global_N_deposition_grid_dataset_2008_2020.rar",
    overwrite = TRUE,
    url = "https://springernature.figshare.com/ndownloader/files/48644623",
    projectTo = to,
    destinationPath = destinationPath,
    fun = "terra::rast"
  )
  # fill the holes
  w <- sum(is.na(values(Ndeposition2)))
  while (w != 0){
    w <- round(sqrt(w))
    # make sure w is a odd number
    if(w %% 2 != 1){
      w = w+1
    }
    Ndeposition2 <- focal(Ndeposition2, w = w, fun = 'mean', na.policy = 'only')
    w <- sum(is.na(values(Ndeposition2)))
  }
  # mask to study area
  Ndeposition2 <- maskTo(Ndeposition2, to) |> round()
  
  # convert from kg/ha/yr to kg/m2/yr
  Ndeposition <- c(Ndeposition1/10000, Ndeposition2/10000) 
  names(Ndeposition) <- as.character(c(year1, year2))
  
  return(Ndeposition)
}

# Get the albedo per land cover
rvestAlbedoTable <- function(destinationPath){
  
  # extract Table 1 of Gao et al., 2005 (https://doi.org/10.1029/2004JD005190)
  url <- "http://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2004JD005190"
  tmpfile <- file.path(destinationPath, "tmp.html")
  download.file(url, destfile = tmpfile, quiet = TRUE)
  GaoEtAl2005 <- read_html(tmpfile)
  file.remove(tmpfile)
  tables <- GaoEtAl2005 %>% html_table(fill = TRUE)
  table1 <- tables[[1]]
  
  # Only keep the snow-free shortwave albedo entries
  shortwaveAlbedo <- table1[c(26:34), c(2,3,5,7,9)]
  
  # Add column names
  names(shortwaveAlbedo) <- c("IGBPclass", "lat6070", "lat5060", "lat4050", "lat3040")
  
  # Fix missing entries, set to NA
  shortwaveAlbedo[shortwaveAlbedo=="â\u0080\u0093"] <- NA
  
  # Convert entries to numeric
  shortwaveAlbedo[, c(2:5)] <- apply(shortwaveAlbedo[, c(2:5)], MARGIN = 2, as.numeric)
  
  return(as.data.frame(shortwaveAlbedo))
}

# Get the co2 concentration for a period for different representative concentration pathways (RCPs)
prepCo2Concentration <- function(firstYear, lastYear, scenario, destinationPath){
  dir.create(file.path(destinationPath, "co2"), showWarnings = FALSE)
  if (!(scenario %in% c("RCP26", "RCP45", "RCP60", "RCP85"))) {
    stop(
      "Available representative concentration pathways for co2 concentration ",
      "are RCP26, RCP45, RCP60, or RCP85"
    )
  }
  targetFile <- paste0(scenario, "_MIDYR_CONC.DAT")
  url <- paste0(
    "https://raw.githubusercontent.com/ClimateChangeEcon/Climate_in_Climate_Economics/refs/heads/main/calibration_data/EmiAndConcData/",
    targetFile
  )
  co2concentrations <- prepInputs(
    targetFile = targetFile,
    url = url,
    fun = {
      co2Data <- read.table(targetFile, skip = 39)
      co2Data <- co2Data[, c(1, 4)]
      co2Data
    },
    destinationPath = destinationPath,
    overwrite = TRUE
  ) |> Cache()
  yearToKeep <- firstYear <= co2concentrations[,1] & co2concentrations[,1] <= lastYear
  co2concentrations <- co2concentrations[yearToKeep, ]
  fileName <- paste("co2", firstYear, lastYear, paste0(scenario, ".txt"), sep = "_")
  CO2write(co2concentrations,
           fileName = file.path(destinationPath, "co2", fileName))
  names(co2concentrations) <- c("year", "co2_ppm")
  return(co2concentrations)
}

#
prepEPC <- function(url, sppEquiv, destinationPath){
  # create a folder for epc in the destinationPath
  dir.create(file.path(destinationPath, "epc"), showWarnings = FALSE)
  # read epc
  epc <- prepInputs(targetFile = "defaultEPC.csv", 
                    url = url, 
                    destinationPath = destinationPath,
                    fun = "data.table::fread",
                    check.names = TRUE,
                    overwrite = TRUE)
  
  # keep only the lines for the species in study
  epc <- epc[epc$speciesId %in% sppEquiv$speciesId, ]
  
  # Write the species-level epcs in the destinationPath/epc folder
  apply(epc, MARGIN = 1, epcWrite2, destinationPath = destinationPath)
  
  return(epc)
}

# Extract the meteorological data
prepClimate <- function(climatePolygons, siteName, simStartYear, simEndYear, nSpinupYears, scenario, climModel, destinationPath){
  # Create a folder where metdata will be saved
  dir.create(file.path(destinationPath, "metdata"), showWarnings = FALSE)
  
  # get latitude and longitude
  latlon <- project(crds(centroids(climatePolygons)), crs(climatePolygons), "EPSG:4326")
  lon <- latlon[,1]
  lat <- latlon[,2]
  
  if ("climatePolygonId" %in% names(climatePolygons)){
    id = climatePolygons$climatePolygonId
  } else {
    id = c(1:length(climatePolygons))
  }
  
  # define the first year to extract climate
  firstYear <- simStartYear - nSpinupYears
  
  # get climate from BioSim
  climate <- generateWeather(
    modelNames = c("Climatic_Daily", "VaporPressureDeficit_Daily"),
    fromYr = firstYear,
    toYr = simEndYear,
    id = id,
    latDeg = lat,
    longDeg = lon,
    rcp = scenario,
    climModel = climModel,
    additionalParms = NULL
  )
  
  # format climate data, 1 ecodistrict at a time
  climOut <- list()
  for (i in id) {
    climate_i <- lapply(climate, FUN = function(x){
      x <- x[x$KeyID == i, ]
      # BiomeBGC do not simulate Feb 29 (always 365 day/yr)
      x <- x[x$Month != 2 | x$Day != 29,]
      x
    })
    
    daylen <- daylength(lat[id == i], 1:365) * 60 * 60
    
    climate_i <- data.frame(
      year = climate_i[["Climatic_Daily"]]$Year,
      yday = 1:365,
      tmax = climate_i[["Climatic_Daily"]]$Tmax,
      tmin = climate_i[["Climatic_Daily"]]$Tmin,
      tday = climate_i[["Climatic_Daily"]]$Tair,
      prcp = climate_i[["Climatic_Daily"]]$Prcp/10, # from mm to cm
      vpd = climate_i[["VaporPressureDeficit_Daily"]]$VaporPressureDeficit * 100, # from hPa to Pa
      srad = climate_i[["Climatic_Daily"]]$SRad,
      daylen = daylen,
      spinup = climate_i[["Climatic_Daily"]]$Year < simStartYear
    )
    
    spinupFileName <- tolower(paste0(
      i,
      "_",
      climModel,
      scenario,
      "_spinup.mtc43"
    ))
    spinupFileName <- file.path(destinationPath, "metdata", spinupFileName)
    
    metWrite(
      metData = climate_i[climate_i$spinup, c(1:10)],
      fileName = spinupFileName,
      siteName  = paste0("Climate Polygon: ", i),
      dataSource = paste(climModel, scenario, sep = ": ")
    )
    
    # Met data for main simulation
    fileName <- tolower(paste0(
      i,
      "_",
      climModel,
      scenario,
      "_",
      simStartYear,
      simEndYear,
      ".mtc43"
    ))
    fileName <- file.path(destinationPath, "metdata", fileName)
    
    metWrite(
      metData = climate_i[!climate_i$spinup, c(1:10)],
      fileName = fileName,
      siteName  = paste0("Climate Polygon: ", i),
      dataSource = paste(climModel, scenario, sep = ": ")
    )
    
    
    climOut[[as.character(i)]] <- climate_i
  }
  
  return(climOut)
}

# Extract elevation raster
prepElevation <- function(studyArea, to){
  # get data from Amazon Web Services Terrain Tiles
  elevation <-  get_elev_raster(
    locations = sf::st_as_sf(studyArea),
    z = 10
  )
  
  # Crop and project to raster
  elevation <- postProcessTo(
    rast(elevation),
    to = to
  )
  
  # round to the nearest 50m
  return(50 * round(elevation/50))
}
