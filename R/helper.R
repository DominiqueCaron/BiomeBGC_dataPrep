### EPC helper function
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

readCO2 <- function(file){
  co2Data <- read.table(file, skip = 39)
  co2Data <- co2Data[ ,c(1,4)]
  return(co2Data)
}

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

### MET helper function
metRead <- function(fileName, nHeaderLines = 4){
  fileName <- "~/repos/BiomeBGCR/inst/inputs/metdata/miss5093.mtc41"
  # Always the same 9 variables
  colNames <- c("year", "yday", "tmax", "tmin", "tday", "prcp", "vpd", "srad", "daylen")
  
  # Read data, skip header
  metData <- read.table(fileName, skip = nHeaderLines, col.names = colNames)
  
  return(metData)
}

metWrite <- function(metData, fileName, studyArea = "XXXX", dataSource = "XXXX"){
  # Create and enter file
  sink(fileName)
  
  # First row with study name and year range
  cat(studyArea, ",", sep = "")
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

## CO2 helper function
CO2Read <- function(fileName){
  fileName <- "~/repos/BiomeBGCR/inst/inputs/co2/co2.txt"
  co2Data <- read.table(fileName, col.names = c("year", "concentration"))
  return(co2Data)
}

CO2write <- function(co2Data, fileName){
  write.table(co2Data, 
              file = fileName, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE)
}

assertEPCproportions <- function(epc){
  # check Litter proportions
  ## If there is only one missing, fill by subtracting the sum of the other two to 1
  epc[, litterProportionN :=
        rowSums(!is.na(.SD)),
      .SDcols = c("litterLabileProportion",
                  "litterCelluloseProportion",
                  "litterLigninProportion")]
  epc[litterProportionN == 2 & is.na(litterLabileProportion), litterLabileProportion := 1 - (litterCelluloseProportion + litterLigninProportion)]
  epc[litterProportionN == 2 & is.na(litterCelluloseProportion), litterCelluloseProportion := 1 - (litterLabileProportion + litterLigninProportion)]
  epc[litterProportionN == 2 & is.na(litterLigninProportion), litterLigninProportion := 1 - (litterCelluloseProportion + litterLabileProportion)]
  ## If the three are provided, make sure that they sum to 1 (scale if necessary)
  epc[, litterProportionSum := litterLabileProportion + litterCelluloseProportion + litterLigninProportion]
  epc[litterProportionN == 3 & litterProportionSum != 1, litterLabileProportion := litterLabileProportion / litterProportionSum]
  epc[litterProportionN == 3 & litterCelluloseProportion != 1, litterCelluloseProportion := litterCelluloseProportion / litterProportionSum]
  epc[litterProportionN == 3 & litterLigninProportion != 1, litterLigninProportion := litterLigninProportion / litterProportionSum]
  epc[ , litterProportionN := NULL]
  epc[ , litterProportionSum := NULL]
  
  # do the same with fine root proportions
  ## If there is only one missing, fill by subtracting the sum of the other two to 1
  epc[, FRProportionN :=
        rowSums(!is.na(.SD)),
      .SDcols = c("fineRootLabileProportion",
                  "fineRootCelluloseProportion",
                  "fineRootLigninProportion")]
  epc[FRProportionN == 2 & is.na(fineRootLabileProportion), fineRootLabileProportion := 1 - (fineRootCelluloseProportion + fineRootLigninProportion)]
  epc[FRProportionN == 2 & is.na(fineRootCelluloseProportion), fineRootCelluloseProportion := 1 - (fineRootLabileProportion + fineRootLigninProportion)]
  epc[FRProportionN == 2 & is.na(fineRootLigninProportion), fineRootLigninProportion := 1 - (fineRootCelluloseProportion + fineRootLabileProportion)]
  ## If the three are provided, make sure that they sum to 1 (scale if necessary)
  epc[, FRProportionSum := fineRootLabileProportion + fineRootCelluloseProportion + fineRootLigninProportion]
  epc[FRProportionN == 3 & FRProportionSum != 1, fineRootLabileProportion := fineRootLabileProportion / FRProportionSum]
  epc[FRProportionN == 3 & FRProportionSum != 1, fineRootCelluloseProportion := fineRootCelluloseProportion / FRProportionSum]
  epc[FRProportionN == 3 & FRProportionSum != 1, fineRootLigninProportion := fineRootLigninProportion / FRProportionSum]
  epc[ , FRProportionN := NULL]
  epc[ , FRProportionSum := NULL]
 
  # and finally with dead wood
  ## If there is only one missing, fill by subtracting the sum of the other two to 1
  epc[, DWProportionN :=
        rowSums(!is.na(.SD)),
      .SDcols = c("deadWoodCelluloseProportion",
                  "deadWoodLigninProportion")]
  epc[DWProportionN == 1 & is.na(deadWoodCelluloseProportion), deadWoodCelluloseProportion := 1 - deadWoodLigninProportion]
  epc[DWProportionN == 1 & is.na(deadWoodLigninProportion), deadWoodLigninProportion := 1 - deadWoodCelluloseProportion]
  ## If the three are provided, make sure that they sum to 1 (scale if necessary)
  epc[, DWProportionSum := deadWoodCelluloseProportion + deadWoodLigninProportion]
  epc[DWProportionN == 2 & DWProportionSum != 1, deadWoodCelluloseProportion := deadWoodCelluloseProportion / DWProportionSum]
  epc[DWProportionN == 2 & DWProportionSum != 1, deadWoodLigninProportion := deadWoodLigninProportion / DWProportionSum]
  epc[ , DWProportionN := NULL]
  epc[ , DWProportionSum := NULL]
  
  return(epc)
}

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
      initialVaportPressureDeficit = numeric(),
      finalVaporPressureDeficit = numeric()
    )
  )
}
