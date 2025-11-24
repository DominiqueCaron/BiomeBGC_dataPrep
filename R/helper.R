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
