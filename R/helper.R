## helper functions:
# epcRead
# epcParseLine
# epcWrite2
# epcWrite
# metRead (not used)
# metWrite
# CO2Read
# CO2Write
# assertEPCproportions
# initiateEPC
# getEPC
# lccToAlbedo

#### EPC helper functions ####
# reads .epc files
epcRead <- function(fileName, readValues = TRUE, readMeta = TRUE){
  # Get path to the default epc files if using one of them.
  if(fileName %in% c("c3grass", "c4grass", "dbf", "dnf", "ebf", "enf", "shrub")){
    fileName <- file.path(
      system.file("inputs", package = "BiomeBGCR"),
      "epc",
      paste0(fileName, ".epc")
    )
  }
  
  if (readValues | readMeta){
    con <- file(fileName, open = "r")
    
    # Read header
    epcHeader <- readLines(con, n = 1)
    
    # Read table
    epcTable <- readLines(con)
    
    # Extract values
    functionalType <- sub("^ECOPHYS\\s{2,}", "", epcHeader)
    epcTableParsed <- sapply(epcTable, FUN = epcParseLine, USE.NAMES = FALSE)
    
    # If asked to read values and metadata columns, return as data.frame with attributes
    if (readValues & readMeta){
      epc <- data.frame(
        description = epcTableParsed[1,],
        unit = epcTableParsed[2,],
        value = as.numeric(epcTableParsed[3,])
      )
      # Add attributes
      attr(epc, "title") <- "ECOPHYS"
      attr(epc, "description") <- functionalType
      
      
    } else if (readValues & !readMeta){
      # If asked to read values only, return a vector of values
      epc <- as.numeric(epcTableParsed[3,])
    } else if (!readValues & readMeta){
      # If asked to read metadata columns only, return a data.frame with 2 columns
      epc <- data.frame(
        description = epcTableParsed[1,],
        unit = epcTableParsed[2,]
      )
    }
    close(con)
  } else {
    # If asked to read nothing, return a warning and an empty vector.
    warning("readValues and readMeta are FALSE: returning an empty vector.")
    epc <- c()
    
  }
  return(epc)
}

# parse a single line in epc files.
epcParseLine <- function(epcLine){
  # The value comes first and is followed by >=2 spaces
  epcValue <- strsplit(epcLine, "\\s{1,}")[[1]][1]
  
  # The unit is inbetween parentheses.
  epcUnit <- strsplit(epcLine, "\\s{1,}")[[1]][2]
  epcUnit <- gsub("\\(|\\)", "", epcUnit)
  
  # The description is everything after
  epcDescription <- strsplit(epcLine, "\\s{1,}")[[1]][-c(1,2)]
  epcDescription <- paste(epcDescription, collapse = " ")
  if(epcUnit == "flag"){
    epcDescription <- gsub(" 0", "; 0", epcDescription)
  }
  
  return(c(epcDescription, epcUnit, epcValue))
}

# writes a .epc file from the ecophysiologicalConstants input
epcWrite2 <- function(epc, destinationPath){
  # Get a template
  epcTemplate <- file.path(
    system.file("inputs", package = "BiomeBGCR"),
    "epc", "enf.epc")
  epcOut <- epcRead(fileName = epcTemplate)
  # set Description and fileName
  attr(epcOut, "description") <- paste0(epc["species"], "-dominated stands")
  fileName <- tolower(paste0(gsub(" ", "", epc["species"]), ".epc"))
  fileName <- file.path(destinationPath, 
                        "epc", 
                        fileName)
  #replace values
  epcOut[, "value"] <- as.numeric(epc[-c(1,2)])
  # write file
  epcWrite(epcOut, fileName = fileName)
}

# writes a .epc file from a formatted data.frame
epcWrite <- function(epc, fileName){
  # Create and enter file
  sink(fileName)
  # Header line
  cat("ECOPHYS", rep("", 6), attr(epc, "description"))
  # Data lines
  for (i in 1:nrow(epc)){
    cat("\n")
    epcValue <- round(epc[i,3], 5)
    npad1 <- 14 - nchar(format(epcValue, scientific = FALSE))
    epcUnit <- epc[i,2]
    epcUnit <- paste0("(",epcUnit,")")
    npad2 <- 10 - nchar(epcUnit)
    epcDescription <- epc[i,1]
    if(grepl("flag", epcUnit)){
      epcDescription <- strsplit(epcDescription, "; ")
      npad3 <- 22 - nchar(epcDescription[[1]][1])
      epcDescription <- paste0(
        epcDescription[[1]][1], 
        paste(rep(" ", npad3), collapse = ""),
        epcDescription[[1]][2])
    }
    if(!grepl("flag", epcUnit) & !grepl("yday", epcUnit)){
      if (epcValue == round(epcValue)){
        epcValue <- paste0(as.character(epcValue), ".0")
        npad1 <- npad1 - 2
      }
    }
    if (grepl("\\*", epcUnit)){
      epcUnit <- gsub("\\(\\*", "\\*\\(", epcUnit)
      npad1 <- npad1 - 1
      npad2 <- npad2 + 1
    }
    cat(format(epcValue, scientific = FALSE), rep("", npad1))
    cat(epcUnit, rep("", npad2))
    cat(epcDescription)
  }
  cat("\n")
  sink()
}

# If there are proportions for all but one component, infer the missing proportion (e.g.,prop3 = 1-(prop1+prop2))
fillMissingEPCProportions <- function(epc, pools = c("fineRoot", "litter", "deadWood")){
  for (pool in pools){
    # get the correct columns
    cols <- grep(paste(pool, "Proportion", sep = ".*"), colnames(epc))
    
    # find the lines with 1 missing value
    nonNAcols <- rowSums(!is.na(epc[,..cols]))
    rowToFill <- which((length(cols) - nonNAcols) == 1)
    
    # fill the filling value with 1 - sum(other values)
    epc[rowToFill, (cols) := {
      x <- .SD
      miss <- is.na(x)
      x[miss] <- 1 - rowSums(x, na.rm = TRUE)
      x
    }, .SDcols = cols]
  }
  return(epc)
}

# If all fractions are available, make sure they sum to 1
scaleEPCProportions <- function(epc, pools = c("fineRoot", "litter", "deadWood")){
  for (pool in pools){
    # get the correct columns
    cols <- grep(paste(pool, "Proportion", sep = ".*"), colnames(epc))
    
    # find the lines that have all values and do not sum to 1
    sumFractions <- rowSums(epc[,..cols])
    rowToScale <- which(!is.na(sumFractions) & sumFractions != 1)
    
    # scale values
    if(length(rowToScale) != 0){
      epc[rowToScale, (cols) := {
        x <- .SD
        x <- x/rowSums(x)
        x
      }, .SDcols = cols]
    }
    
  }
  return(epc)
}

# Makes sure that the sums of proportions in epc sum to 1
cleanEPC <- function(epc){
  # scale proportions to 1
  epc <- scaleEPCProportions(epc)
  # round constants
  epc[, "leafAndFineRootTurnover"] <- round(epc[, "leafAndFineRootTurnover"], 3)
  allocationCols <- grep("new", names(epc))
  epc[, allocationCols] <- round(epc[, ..allocationCols], 2)
  CtoNCols <- grep("CtoN", names(epc))
  epc[, CtoNCols] <- round(epc[, ..CtoNCols], 2)
  proportionCols <- grep("Proportion", names(epc))[-1]
  epc[, proportionCols] <- round(epc[, ..proportionCols], 2)
  epc[, c(32:39)] <- round(epc[, c(32:39)], 3)
  
  # After rounding, make sure proportions equal to 1
  epc[, litterLigninProportion := 1-(litterLabileProportion+litterCelluloseProportion)]
  epc[, fineRootLigninProportion  := 1-(fineRootCelluloseProportion +fineRootLabileProportion)]
  epc[, deadWoodLigninProportion := 1-deadWoodCelluloseProportion]
  return(epc)
}

# Initiates an epc data.frame with correct column names
initiateEPC <- function() {
  return(
    data.frame(
      taxa = character(),
      level = factor(levels = c("species", "genus", "pft")),
      woody = integer(),
      evergreen = integer(),
      c3psn = integer(),
      modelPhenology = integer(),
      yearDayToStartGrowth = integer(),
      yearDayToEndLitterfall = integer(),
      transferGrowthPeriod = numeric(),
      litterfallPeriod = numeric(),
      leafAndFineRootTurnover = numeric(),
      annualLiveWoodTurnover = numeric(),
      wholePlantMortality = numeric(),
      fireMortality = numeric(),
      newFineRootCToNewLeafC = numeric(),
      newStemCToNewLeafC = numeric(),
      newLiveWoodCToNewTotalWoodC = numeric(),
      newCoarseRootCToNewStemC = numeric(),
      currentGrowthProportion = numeric(),
      CtoNLeaves = numeric(),
      CtoNLitter = numeric(),
      CtoNFineRoots = numeric(),
      CtoNLiveWood = numeric(),
      CtoNDeadWood = numeric(),
      litterLabileProportion = numeric(),
      litterCelluloseProportion = numeric(),
      litterLigninProportion = numeric(),
      fineRootLabileProportion = numeric(),
      fineRootCelluloseProportion = numeric(),
      fineRootLigninProportion = numeric(),
      deadWoodCelluloseProportion = numeric(),
      deadWoodLigninProportion = numeric(),
      canopyWaterInterception = numeric(),
      canopyLightExtinction = numeric(),
      allToProjectLAI = numeric(),
      canopyAvgSLA = numeric(),
      shadedSLAToSunlightSLA = numeric(),
      NfractionInRuBisCO = numeric(),
      maxStomatalConductance = numeric(),
      cuticularConductance = numeric(),
      boundaryLayerConductance = numeric(),
      initialLeafWaterPotential = numeric(),
      finalLeafWaterPotential = numeric(),
      initialVaporPressureDeficit = numeric(),
      finalVaporPressureDeficit = numeric()
    )
  )
}

# Get the ecophysiological constants for species in sppEquiv. Starts by retrieving
# species-level epc, then tries genus-level, and finally plant functional type level.
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

#### MET helper functions ####

# Reads .mtc43 data
# Currently not used.
metRead <- function(fileName, nHeaderLines = 4){
  # Always the same 9 variables
  colNames <- c("year", "yday", "tmax", "tmin", "tday", "prcp", "vpd", "srad", "daylen")
  
  # Read data, skip header
  metData <- read.table(fileName, skip = nHeaderLines, col.names = colNames)
  
  return(metData)
}

# Write .mtc43 data
metWrite <- function(metData, fileName, siteName = "XXXX", dataSource = "XXXX"){
  # Create and enter file
  sink(fileName)
  
  # First row with study name and year range
  cat(siteName, ",", sep = "")
  cat(paste(range(metData$year), collapse = "-"))
  cat("\n")
  
  # Second row on information about the source of the data
  cat(dataSource, "OUTPUT FILE :", format(Sys.time(), "%a %b %d %X %Y"))
  cat("\n")
  
  # Third row is the variable names (always the same)
  cat("  year  yday    Tmax    Tmin    Tday    prcp      VPD     srad  daylen")
  cat("\n")
  
  # Fourth row is the variable units (always the same)
  cat("             (deg C) (deg C) (deg C)    (cm)     (Pa)  (W m-2)     (s)")
  cat("\n")
  
  # The rest are the data.
  # Following a template to make a clean file
  for (i in 1:nrow(metData)){
    # 6 characters for year and yday
    cat(formatC(metData[i, 1], width = 6))
    cat(formatC(metData[i, 2], width = 6))
    # 8 character for Tmax, Tmin, Tday, and prcp
    cat(formatC(metData[i, 3], digits = 2, format = "f", width = 8))
    cat(formatC(metData[i, 4], digits = 2, format = "f", width = 8))
    cat(formatC(metData[i, 5], digits = 2, format = "f", width = 8))
    cat(formatC(metData[i, 6], digits = 2, format = "f", width = 8))
    # 9 characters for prcp and VPD
    cat(formatC(metData[i, 7], digits = 2, format = "f", width = 9))
    cat(formatC(metData[i, 8], digits = 2, format = "f", width = 9))
    # 8 characters for daylen
    cat(formatC(metData[i, 9], format = "d", width = 8))
    cat("\n")
  }
  sink()
}

#### CO2 helper function ####

# reads co2 concentration data file
# currently not used
CO2Read <- function(fileName){
  fileName <- "~/repos/BiomeBGCR/inst/inputs/co2/co2.txt"
  co2Data <- read.table(fileName, col.names = c("year", "concentration"))
  return(co2Data)
}

# writes co2 concentration data file
CO2write <- function(co2Data, fileName){
  write.table(co2Data, 
              file = fileName, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE)
}

#### Albedo helper functions ####

# Gets the soil albedo given the land cover class and latitude. 
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

#### Biome-BGC output helper functions ####
getOutputDescription <- function(outputId){
  # get the output and their id from github
  url <- "https://raw.githubusercontent.com/bpbond/Biome-BGC/refs/heads/master/src/bgclib/output_map_init.c"
  txt <- readLines(url)
  
  # extract the relevant lines
  name_lines <- txt[grepl("output_map\\[[0-9]", txt)]
  
  outputTbl <- data.table(
    output_id = as.integer(sub(".*\\[([^]]+)\\].*", "\\1", name_lines)),
    output_name = sub(".*&([^->]+)->([^;]+);.*", "\\1.\\2", name_lines)
  )
  
  out <- outputTbl[output_id %in% outputId]$output_name
  
  return(out)
}
