## Functions used to in the .inputObjects event:
# prepNTEMSDominantSpecies
# prepSoilTexture
# prepNdeposition
# rvestAlbedoTable
# prepCo2Concentration
# prepWhite2000Table
# prepWhite2000EPC
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
    method = method
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
prepSoilTexture <- function(destinationPath, to){
  sand <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Sand/Sand_X15_30_cm_100m1980-2000v1.tif",
    targetFile = "Sand_X15_30_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to
  ) |> round(digits = -1)
  clay <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Clay/Clay_X15_30_cm_100m1980-2000v1.tif",
    targetFile = "Clay_X15_30_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    to = to
  ) |> round(digits = -1)
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
    destinationPath = destinationPath,
    to = to,
    fun = "terra::rast"
  ) |> round()
  Ndeposition2 <- prepInputs(
    targetFile = paste0("mean_totN_", year2, "_hm.tif"),
    archive = "Global_N_deposition_grid_dataset_2008_2020.rar",
    overwrite = TRUE,
    url = "https://springernature.figshare.com/ndownloader/files/48644623",
    destinationPath = destinationPath,
    to = to,
    fun = "terra::rast"
  ) |> round()
  Ndeposition <- c(Ndeposition1/10000, Ndeposition2/10000) # convert from kg/ha/yr to kg/m2/yr
  names(Ndeposition) <- as.character(c(year1, year2))
  
  return(Ndeposition)
}

# Get the albedo per land cover
rvestAlbedoTable <- function(){
  
  # extract Table 1 of Gao et al., 2005 (https://doi.org/10.1029/2004JD005190)
  GaoEtAl2005 <- read_html("http://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2004JD005190")
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
    destinationPath = destinationPath
  ) |> Cache()
  yearToKeep <- firstYear <= co2concentrations[,1] & co2concentrations[,1] <= lastYear
  co2concentrations <- co2concentrations[yearToKeep, ]
  fileName <- paste("co2", firstYear, lastYear, paste0(scenario, ".txt"), sep = "_")
  CO2write(co2concentrations,
           fileName = file.path(destinationPath, "co2", fileName))
  names(co2concentrations) <- c("year", "co2_ppm")
  return(co2concentrations)
}

# Clean and format tables from White et al., 2000 (EPC constant)
prepWhite2000Table <- function(tbl, value.var){
  # Remove plant functional types
  tbl <- tbl[tbl$Foliage.Nature %in% c("Evergreen needle leaf forest", "Deciduous broad leaf forest"),]
  # Remove extra-spaces
  tbl$Species <- trimws(tbl$Species)
  # Standardize genus-level taxa
  tbl$Species <- gsub(pattern = " spp.", "", tbl$Species)
  # Transform into wide format
  tbl <- dcast(as.data.table(tbl), Species ~ Parameter, value.var = value.var, fun.aggregate = mean)
  # Extract genus
  tbl[, Genus := sub(" .*", "", Species)]
  # Table of species-level constants
  tbl_sp <- tbl[Species != Genus, ]
  tbl_sp[, Genus := NULL]
  tbl_sp[, Level := "Species"]
  # Table of genus-level constants - take mean across species
  tbl[, Species := NULL]
  tbl_genus <- tbl[, lapply(.SD, mean, na.rm = TRUE), by = Genus]
  tbl_genus[, Level := "Genus"]
  setnames(tbl_genus, "Genus", "Species")
  # Combine genus and species-level constants
  out <- rbindlist(list(tbl_sp, tbl_genus))
  return(out)
}

# Prepares epc from the tables in white et al., 2000
prepWhite2000EPC <- function(url, sppEquiv, destinationPath){
  # create a folder for epc in the destinationPath
  dir.create(file.path(destinationPath, "epc"), showWarnings = FALSE)
  
  # get a epc template file
  epc <- initiateEPC()
  
  ## append the enf and dbf plant functional types
  epc <- rbindlist(list(epc, data.table(taxa = c("enf", "dbf"), level = "pft")), fill = TRUE)
  enf <- epcRead("enf")
  dbf <- epcRead("dbf")
  epc[1, c(3:45)] <- data.table(matrix(enf$value, 1, 43))
  epc[2, c(3:45)] <- data.table(matrix(dbf$value, 1, 43))
  
  ## Prepare parameters from White 2000
  # Table 1
  white2000_1 <- prepInputs(targetFile = "White_spreadsheet1.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  # A bit of cleaning, removing typos and lines we do not need.
  white2000_1$Species[white2000_1$Species == "Prunus pennsylvanica"] <- "Prunus pensylvanica"
  white2000_1$Species[white2000_1$Species == "Rubus alleghaniensis"] <- "Rubus allegheniensis"
  white2000_1$Species[white2000_1$Species == "Rubus allighaniensis"] <- "Rubus allegheniensis"
  white2000_1$Species[white2000_1$Species == "Betula Papyrifera"] <- "Betula papyrifera"
  white2000_1$Species[white2000_1$Species == "Sequoiadendron gigant."] <- "Sequoiadendron giganteum"
  white2000_1$Species[white2000_1$Species == "Tillia americana"] <- "Tilia americana"
  white2000_1 <- white2000_1[!(white2000_1$Species %in% c("DBF", "ENF", "Lonicera x bella", "Mixed deciduous", "Mixed pine", "Nyssa-Acer", "Quercus prinus/rubra", "Rain forest", "Tropical deciduous forest", "Tsuga/Picea", "Wood", "Picea/Abies", "Cedar")), ]
  white2000_1 <- prepWhite2000Table(white2000_1, value.var = "Value")
  names(white2000_1) <- c(
    "taxa",
    "CtoNDeadWood",
    "CtoNFineRoots",
    "newLiveWoodCToNewTotalWoodC",
    "allToProjectLAI",
    "CtoNLeaves",
    "leafAndFineRootTurnover",
    "canopyLightExtinction",
    "CtoNLitter",
    "newCoarseRootCToNewStemC",
    "newFineRootCToNewLeafC",
    "newStemCToNewLeafC",
    "canopyAvgSLA",
    "level"
  )
  
  # White 2000 table 2
  white2000_2 <- prepInputs(targetFile = "White_spreadsheet2.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  # a bit of cleaning, removing typos and lines we do not need.
  white2000_2$Species[white2000_2$Species == "Abies concolr"] <- "Abies concolor"
  white2000_2$Species[white2000_2$Species == "Prunus pensylvannica"] <- "Prunus pensylvanica"
  white2000_2$Species <- gsub("Abiea", "Abies",  white2000_2$Species)
  white2000_2 <- white2000_2[!(white2000_2$Species == "ENF"), ]
  white2000_2[white2000_2 == -999] <- NA
  white2000_2[, c("Labile", "Cellulose", "Lignin")] <- white2000_2[, c("Labile", "Cellulose", "Lignin")]/100
  white2000_2 <- prepWhite2000Table(white2000_2, value.var = c("Labile", "Cellulose", "Lignin"))
  names(white2000_2) <- c(
    "taxa",
    "fineRootLabileProportion",
    "litterLabileProportion",
    "fineRootCelluloseProportion",
    "litterCelluloseProportion",
    "fineRootLigninProportion",
    "litterLigninProportion",
    "level")
  # If there are proportions for all but one component, infer the missing proportion (e.g.,prop3 = 1-(prop1+prop2))
  white2000_2 <- fillMissingEPCProportions(white2000_2, pools = c("fineRoot", "litter"))
  white2000 <- merge(white2000_1, white2000_2, by = c("taxa", "level"), all = TRUE)
  
  # White 2000 table 3
  white2000_3 <- prepInputs(targetFile = "White_spreadsheet3.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  white2000_3$Species[white2000_3$Species == "Pinus elliotii"] <- "Pinus elliottii"
  white2000_3$Species[white2000_3$Species == "Pinus Taeda"] <- "Pinus taeda"
  white2000_3$Species[white2000_3$Species == "Robinea pseudoacacia"] <- "Robinia pseudoacacia"
  white2000_3[white2000_3 == -999] <- NA # These are incorrect
  # except for the lines with NAs, the column are mixed-up
  white2000_3[!is.na(white2000_3$Cellulose), c("Lignin", "Cellulose")] <-  white2000_3[!is.na(white2000_3$Cellulose), c("Cellulose", "Lignin")]
  white2000_3[, c("Lignin", "Cellulose")] <- white2000_3[, c("Lignin", "Cellulose")]/100 
  white2000_3 <- prepWhite2000Table(white2000_3, value.var = c("Lignin", "Cellulose"))
  names(white2000_3) <- c(
    "taxa",
    "deadWoodLigninProportion",
    "deadWoodCelluloseProportion",
    "level")
  white2000_3 <- fillMissingEPCProportions(white2000_3, pools = "deadWood")
  white2000 <- merge(white2000, white2000_3, by = c("taxa", "level"), all = TRUE)
  
  # White 2000 table 4
  white2000_4 <- prepInputs(targetFile = "White_spreadsheet4.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  white2000_4[white2000_4 == -999] <- NA
  white2000_4 <- prepWhite2000Table(white2000_4, value.var = c("Initial", "Final"))
  names(white2000_4) <- c(
    "taxa",
    "initialLeafWaterPotential",
    "finalLeafWaterPotential",
    "level"
  )
  white2000 <- merge(white2000, white2000_4, by = c("taxa", "level"), all = TRUE)
  
  # White 2000 table 5
  white2000_5 <- prepInputs(targetFile = "White_spreadsheet5.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  white2000_5[,c("Initial", "Final")] <- white2000_5[,c("Initial", "Final")] * 1000 # from kPa to Pa
  white2000_5 <- prepWhite2000Table(white2000_5, value.var = c("Initial", "Final"))
  names(white2000_5) <- c(
    "taxa",
    "initialVaporPressureDeficit",
    "finalVaporPressureDeficit",
    "level"
  )
  white2000 <- merge(white2000, white2000_5, by = c("taxa", "level"), all = TRUE)
  
  epc <- rbindlist(list(epc, white2000), fill = TRUE)
  # check that the sums of proportions equal to 1
  epc <- scaleEPCProportions(epc)
  # Get the correct epc for the species in sppEquiv
  epc <- getEPC(epc, sppEquiv)
  epc <- cleanEPC(epc)
  # Write the species-level epcs in the destinationPath/epc folder
  apply(epc, MARGIN = 1, epcWrite2, destinationPath = destinationPath)
  return(epc)
}

# Extract the meteorological data
prepClimate <- function(climatePolygons, siteName, firstYear, lastYear, lastSpinupYear, scenario, climModel, destinationPath){
  # Create a folder where metdata will be saved
  dir.create(file.path(destinationPath, "metdata"), showWarnings = FALSE)
  
  # get latitude and longitude
  latlon <- project(crds(centroids(climatePolygons)), crs(climatePolygons), "EPSG:4326")
  lon <- latlon[,1]
  lat <- latlon[,2]
  
  if ("ECODISTRIC" %in% names(climatePolygons)){
    id = climatePolygons$ECODISTRIC
  } else {
    id = c(1:length(climatePolygons))
  }
  
  # get climate from BioSim
  climate <- generateWeather(
    modelNames = c("Climatic_Daily", "VaporPressureDeficit_Daily"),
    fromYr = firstYear,
    toYr = lastYear,
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
      spinup = climate_i[["Climatic_Daily"]]$Year <= lastSpinupYear
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
      firstYear,
      lastYear,
      ".mtc43"
    ))
    fileName <- file.path(destinationPath, "metdata", fileName)
    
    metWrite(
      metData = climate_i[, c(1:10)],
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
  
  # round to 10m
  return(round(elevation, -1))
}
