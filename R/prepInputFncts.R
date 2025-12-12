prepEPC <- function(dataSource, destinationPath, to = NULL){
  dir.create(file.path(destinationPath, "epc"), showWarnings = FALSE)
  # Option 1: Getting epc from default functional types and based on landcover. 
  if (dataSource == "functionalTypes"){
    lcc <- prepInputs(url = "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/SCANFI/v1/SCANFI_att_nfiLandCover_SW_2020_v1.2.tif",
                      destinationPath= destinationPath,
                      cropTo = buffer(to, 30),
                      projectTo = crs(to)) |> Cache()
    lcc <- extract(lcc, to)
    functionalTypes <- lccToFuncType(lcc[, 2])
    epc <- lapply(functionalTypes, prepEPC, destinationPath)
    
  } else {
    # Option 2: User sets the functional type to uses
    if(dataSource %in% c("c3grass", "c4grass", "dbf", "dnf", "ebf", "enf", "shrub")){
      fileName <- file.path(
        system.file("inputs", package = "BiomeBGCR"),
        "epc",
        paste0(dataSource, ".epc")
      )
    } else {
      # Option 3: User provides the path to the epc file.
      fileName <- dataSource
    }
    # Copy the epc file in the input folder
    newFileName <- file.path(destinationPath, "epc", basename(fileName))
    file.copy(fileName, 
              newFileName,
              overwrite = T
    )
    epc <- epcRead(newFileName)
  }
  return(epc)
}


prepSoilTexture <- function(destinationPath, studyArea){
  sand <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Sand/Sand_X15_30_cm_100m1980-2000v1.tif",
    targetFile = "Sand_X15_30_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    cropTo = buffer(studyArea, 100),
    projectTo = crs(studyArea)
  )
  clay <- prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/psm/Clay/Clay_X15_30_cm_100m1980-2000v1.tif",
    targetFile = "Clay_X15_30_cm_100m1980-2000v1.tif",
    destinationPath = destinationPath,
    cropTo = buffer(studyArea, 100),
    projectTo = crs(studyArea)
  )
  silt <- 100 - (sand + clay)
  soilTexture <- c(sand, silt, clay)
  names(soilTexture) <- c("sand", "silt", "clay")
  return(soilTexture)
}

prepNdeposition <- function(destinationPath, studyArea, year1, year2){
  Ndeposition1 <- prepInputs(
    targetFile = paste0("mean_totN_", year1, "_hm.tif"),
    archive = "Global_N_deposition_grid_dataset_2008_2020.rar",
    overwrite = TRUE,
    url = "https://springernature.figshare.com/ndownloader/files/48644623",
    destinationPath = destinationPath,
    cropTo = buffer(studyArea, 1000),
    projectTo = crs(studyArea),
    fun = "terra::rast"
  )
  Ndeposition2 <- prepInputs(
    targetFile = paste0("mean_totN_", year2, "_hm.tif"),
    archive = "Global_N_deposition_grid_dataset_2008_2020.rar",
    overwrite = TRUE,
    url = "https://springernature.figshare.com/ndownloader/files/48644623",
    destinationPath = destinationPath,
    cropTo = buffer(studyArea, 1000),
    projectTo = crs(studyArea),
    fun = "terra::rast"
  )
  Ndeposition <- c(Ndeposition1/10000, Ndeposition2/10000) # convert from kg/ha/yr to kg/m2/yr
  names(Ndeposition) <- as.character(c(year1, year2))
  return(Ndeposition)
}

# NFI land cover class values:
# 1 = Bryoid
# 2 = Herbs
# 3 = Rock
# 4 = Shrub
# 5 = Tree cover is mainly “Treed broadleaf”
# 6 = Tree cover is mainly “Treed conifer”
# 7 = Tree cover is mainly “Treed mixed”
# 8 = Water
lccToFuncType <- function(lcc){
  if (any(lcc %in% c(1,3,8))){
    stop("One of the study site is bryoids, rock or water.")
  }
  # All herbs are set as c3grass- TODO: is that ok??
  funcType <- rep("c3grass", length(lcc))
  funcType[lcc == 4] <- "shrub"
  funcType[lcc == 5] <- "dbf"
  funcType[lcc == 6] <- "enf"
  # Mixed forests?? set as evergreen needleleaf forest. TODO: is that ok??
  funcType[lcc == 7] <- "enf"
}

# NFI land cover class values:
# 1 = Bryoid
# 2 = Herbs
# 3 = Rock
# 4 = Shrub
# 5 = Tree cover is mainly “Treed broadleaf”
# 6 = Tree cover is mainly “Treed conifer”
# 7 = Tree cover is mainly “Treed mixed”
# 8 = Water
lccToAlbedo <- function(lcc, albedoTable, studyArea){
  if (any(lcc %in% c(1,3,8))){
    stop("One of the study site is bryoids, rock or water.")
  }
  latitude <- crds(project(studyArea, "+proj=longlat +ellps=WGS84 +datum=WGS84"))[2]
  if(latitude > 60){
    colId <- 3
  } else if (latitude > 50){
    colId <- 4
  } else if (latitude > 40){
    colId <- 5
  } else {
    colId <- 6
  }
  albedo <- rep(0.2, length(lcc))
  albedo[lcc == 2] <- albedoTable[7, colId]
  albedo[lcc == 4] <- albedoTable[5, colId]
  albedo[lcc == 5] <- albedoTable[3, colId]
  albedo[lcc == 6] <- albedoTable[1, colId]
  albedo[lcc == 7] <- albedoTable[4, colId]
  return(round(albedo, digits = 2))
}

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
  co2concentration <- prepInputs(targetFile = targetFile,
                                 url = url,
                                 fun = readCO2,
                                 destinationPath = destinationPath)
  yearToKeep <- firstYear <= co2concentrations[,1] & co2concentrations[,1] <= lastYear
  co2concentrations <- co2concentrations[yearToKeep, ]
  fileName <- paste("co2", firstYear, lastYear, paste0(scenario, ".txt"), sep = "_")
  write.table(co2concentrations, file = file.path(destinationPath, "co2", fileName))
  names(co2concentrations) <- c("year", "co2_ppm")
  return(co2concentrations)
}

prepWhite2010Table <- function(tbl, value.var){
  tbl <- tbl[tbl$Foliage.Nature %in% c("Evergreen needle leaf forest", "Deciduous broad leaf forest"),]
  tbl$Species <- trimws(tbl$Species)
  tbl$Species <- gsub(pattern = " spp.", "", tbl$Species)
  tbl <- dcast(as.data.table(tbl), Species ~ Parameter, value.var = value.var, fun.aggregate = mean)
  tbl[, Genus := sub(" .*", "", Species)]
  tbl_sp <- tbl[Species != Genus, ]
  tbl_sp[, Genus := NULL]
  tbl_sp[, Level := "Species"]
  tbl[, Species := NULL]
  tbl_genus <- tbl[, lapply(.SD, mean, na.rm = TRUE), by = Genus]
  tbl_genus[, Level := "Genus"]
  setnames(tbl_genus, "Genus", "Species")
  out <- rbindlist(list(tbl_sp, tbl_genus))
  return(out)
}

prepWhite2010EPC <- function(url, sppEquiv, destinationPath){
  dir.create(file.path(destinationPath, "epc"), showWarnings = FALSE)
  epc <- initiateEPC()
  
  ## append the enf and dbf plant functional types
  epc <- rbindlist(list(epc, data.table(taxa = c("enf", "dbf"), level = "pft")), fill = TRUE)
  enf <- epcRead("enf")
  dbf <- epcRead("dbf")
  epc[1, c(3:45)] <- data.table(matrix(enf$value, 1, 43))
  epc[2, c(3:45)] <- data.table(matrix(dbf$value, 1, 43))
  
  ## Prepare parameters from White 2010
  # Table 1
  white2010_1 <- prepInputs(targetFile = "White_spreadsheet1.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  # A bit of cleaning, removing typos and lines we do not need.
  white2010_1$Species[white2010_1$Species == "Prunus pennsylvanica"] <- "Prunus pensylvanica"
  white2010_1$Species[white2010_1$Species == "Rubus alleghaniensis"] <- "Rubus allegheniensis"
  white2010_1$Species[white2010_1$Species == "Rubus allighaniensis"] <- "Rubus allegheniensis"
  white2010_1$Species[white2010_1$Species == "Betula Papyrifera"] <- "Betula papyrifera"
  white2010_1$Species[white2010_1$Species == "Sequoiadendron gigant."] <- "Sequoiadendron giganteum"
  white2010_1$Species[white2010_1$Species == "Tillia americana"] <- "Tilia americana"
  white2010_1 <- white2010_1[!(white2010_1$Species %in% c("DBF", "ENF", "Lonicera x bella", "Mixed deciduous", "Mixed pine", "Nyssa-Acer", "Quercus prinus/rubra", "Rain forest", "Tropical deciduous forest", "Tsuga/Picea", "Wood", "Picea/Abies", "Cedar")), ]
  white2010_1 <- prepWhite2010Table(white2010_1, value.var = "Value")
  names(white2010_1) <- c(
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
  
  # White 2010 table 2
  white2010_2 <- prepInputs(targetFile = "White_spreadsheet2.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  # a bit of cleaning, removing typos and lines we do not need.
  white2010_2$Species[white2010_2$Species == "Abies concolr"] <- "Abies concolor"
  white2010_2$Species[white2010_2$Species == "Prunus pensylvannica"] <- "Prunus pensylvanica"
  white2010_2$Species <- gsub("Abiea", "Abies",  white2010_2$Species)
  white2010_2 <- white2010_2[!(white2010_2$Species == "ENF"), ]
  white2010_2[white2010_2 == -999] <- NA
  white2010_2[, c("Labile", "Cellulose", "Lignin")] <- white2010_2[, c("Labile", "Cellulose", "Lignin")]/100
  white2010_2 <- prepWhite2010Table(white2010_2, value.var = c("Labile", "Cellulose", "Lignin"))
  names(white2010_2) <- c(
    "taxa",
    "fineRootLabileProportion",
    "litterLabileProportion",
    "fineRootCelluloseProportion",
    "litterCelluloseProportion",
    "fineRootLigninProportion",
    "litterLigninProportion",
    "level")
  white2010 <- merge(white2010_1, white2010_2, by = c("taxa", "level"), all = TRUE)
  
  # White 2010 table 3
  white2010_3 <- prepInputs(targetFile = "White_spreadsheet3.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  white2010_3$Species[white2010_3$Species == "Pinus elliotii"] <- "Pinus elliottii"
  white2010_3$Species[white2010_3$Species == "Pinus Taeda"] <- "Pinus taeda"
  white2010_3[white2010_3 == -999] <- NA
  white2010_3[, c("Lignin", "Cellulose")] <- white2010_3[, c("Lignin", "Cellulose")]/100
  white2010_3 <- prepWhite2010Table(white2010_3, value.var = c("Lignin", "Cellulose"))
  names(white2010_3) <- c(
    "taxa",
    "deadWoodLigninProportion",
    "deadWoodCelluloseProportion",
    "level")
  white2010 <- merge(white2010, white2010_3, by = c("taxa", "level"), all = TRUE)
  
  # White 2010 table 4
  white2010_4 <- prepInputs(targetFile = "White_spreadsheet4.csv", 
                            url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link", 
                            destinationPath = destinationPath,
                            fun = "data.table::fread",
                            check.names = TRUE)
  white2010_4[white2010_4 == -999] <- NA
  white2010_4 <- prepWhite2010Table(white2010_4, value.var = c("Initial", "Final"))
  names(white2010_4) <- c(
    "taxa",
    "initialLeafWaterPotential",
    "finalLeafWaterPotential",
    "level"
  )
  white2010 <- merge(white2010, white2010_4, by = c("taxa", "level"), all = TRUE)
  
  epc <- rbindlist(list(epc, white2010), fill = TRUE)
  epc <- assertEPCproportions(epc)
  epc <- getEPC(epc, sppEquiv)
  apply(epc, MARGIN = 1, epcWrite2, destinationPath = destinationPath)
  return(epc)
}

prepClimate <- function(studyArea, siteName, firstYear, lastYear, scenario, climModel){
  # get latitude and longitude
  latlon <- round(crds(project(studyArea, "+proj=longlat +ellps=WGS84 +datum=WGS84")), 2)
  lon <- latlon[1]
  lat <- latlon[2]
  
  # get climate
  climate <- generateWeather(c("Climatic_Daily", "VaporPressureDeficit_Daily"), firstYear, lastYear, siteName, lat, lon, rcp = scenario, climModel = climModel, additionalParms = NULL)
  
  # fix leap years and organize columns
  
  # get day length
  daylen <- daylength(lat, 1:365)

    
}

getEPC <- function(epcTable, sppEquiv){
  epcSpecies <- epcTable[level == "Species"]
  epcGenus <- epcTable[level == "Genus", ]
  epcPFT <- epcTable[level == "pft", ]
  epc_cols <- colnames(epcTable[,-c(1,2)])
  res <- as.data.table(sppEquiv)
  
  # 1. Species-level traits
  res <- merge(res, epcSpecies, by.x = "species", by.y = "taxa", all.x = TRUE)
  
  # 2. Genus-level (adds trait.genus)
  res <- merge(res, epcGenus, by.x = "genus", by.y = "taxa", all.x = TRUE, suffixes = c("", "_genus"))
  
  # 3. PFT-level (adds trait.pft)
  res <- merge(res, epcPFT, by.x = "PFT", by.y = "taxa", all.x = TRUE, suffixes = c("", "_pft"))
  
  for (col in epc_cols) {
    res[, (col) := fifelse(!is.na(get(col)), 
                           get(col), 
                           fifelse(!is.na(get(paste0(col, "_genus"))),
                                   get(paste0(col, "_genus")),
                                   get(paste0(col, "_pft"))
                           )
    )
    ]
  }
  cols <- c("speciesId", "species", epc_cols)
  res <- res[ , ..cols]
  return(res)
}