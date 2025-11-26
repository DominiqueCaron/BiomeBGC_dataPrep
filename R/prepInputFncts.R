prepEPC <- function(dataSource, destinationPath, to = NULL){
  dir.create(file.path(destinationPath, "epc"), showWarnings = FALSE)
  # Option 1: Getting epc from default functional types and based on landcover. 
  if (dataSource == "functionalTypes"){
    lcc <- prepInputs(url = "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/SCANFI/v1/SCANFI_att_nfiLandCover_SW_2020_v1.2.tif",
                      destinationPath= destinationPath,
                      cropTo = buffer(to, 30),
                      projectTo = crs(to))
    lcc <- extract(lcc, to)
    functionalTypes <- lccToFuncType(lcc[,2])
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
