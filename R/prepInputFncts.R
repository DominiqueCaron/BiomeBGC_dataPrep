prepEPC <- function(dataSource, destinationPath, to = NULL){
  dir.create(file.path(destinationPath, "epc"), showWarnings = FALSE)
  # Option 1: Getting epc from default functional types and based on landcover. 
  if (dataSource == "functionalTypes"){
    lcc <- prepInputs(url = "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/SCANFI/v1/SCANFI_att_nfiLandCover_SW_2020_v1.2.tif",
                      destinationPath= destinationPath,
                      cropTo = buffer(to, 30),
                      projectTo = crs(to))
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

prepNFixationRate <- function()

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

getNfixationRate <- function(studyArea, destinationPath){
  # Get data from https://doi.org/10.5066/P13THKNR
  # There is a high-res and low-res possibility
  # Using low-res for now
  NfixationRates <- prepInputs(
    targetFile = "BNF_total_central_1.tif",
    overwrite = TRUE,
    url = "https://www.sciencebase.gov/catalog/file/get/66b53cc6d34eebcf8bb3850e?f=__disk__67%2Fdf%2F6a%2F67df6a59f896d547205ddb20da99ec72db7a6b10",
    destinationPath = destinationPath,
    cropTo = buffer(studyArea, 1000),
    projectTo = crs(studyArea),
    fun = "terra::rast"
  )
  NfixationRates <- extract(NfixationRates, studyArea)
  NfixationRates <- NfixationRates[,2] / 10000 # convert from kg/ha/yr to kg/m2/yr
  return(round(NfixationRates, 4))
}

getNdeposition <- function(studyArea, year, destinationPath){
  # Get data from https://www.nature.com/articles/s41467-024-55606-y
  # Available years: 2008 to 2020
  Ndeposition <- prepInputs(
    targetFile = paste0("mean_totN_", year, "_hm.tif"),
    archive = "Global_N_deposition_grid_dataset_2008_2020.rar",
    overwrite = TRUE,
    url = "https://springernature.figshare.com/ndownloader/files/48644623",
    destinationPath = destinationPath,
    cropTo = buffer(studyArea, 1000),
    projectTo = crs(studyArea),
    fun = "terra::rast"
  )
  Ndeposition <- extract(Ndeposition, studyArea)
  Ndeposition <- Ndeposition[,2] / 10000 # convert from kg/ha/yr to kg/m2/yr
  return(round(Ndeposition, 4))
}

getAlbedo <- function(studyArea, year, destinationPath){

}
