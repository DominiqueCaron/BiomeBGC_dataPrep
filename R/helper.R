epcRead <- function(fileName){
  con <- file(fileName, open = "r")
  
  # Read header
  epcHeader <- readLines(con, n = 1)
  
  # Read table
  epcTable <- readLines(con)
  
  # Extract values
  functionalType <- sub("^ECOPHYS\\s{2,}", "", epcHeader)
  epcTableParsed <- sapply(epcTable, FUN = epcParseLine, USE.NAMES = FALSE)
  
  # Get into data.frame
  epcDF <- data.frame(
    description = epcTableParsed[1,],
    unit = epcTableParsed[2,],
    value = as.numeric(epcTableParsed[3,])
  )
  
  # Add attributes
  attr(epcDF, "title") <- "ECOPHYS"
  attr(epcDF, "description") <- functionalType
  
  return(epcDF)
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
  fileName = "test.epc"
  epc <- c3grass
  # Create and enter file
  sink(fileName)
  # Header line
  cat("ECOPHYS", rep("", 6), attr(epc, "description"))
  # Data lines
  for (i in 1:nrow(epc)){
    cat("\n")
    epcValue <- round(epc[i,3], 5)
    npad1 <- 14 - nchar(as.character(epcValue))
    epcUnit <- epc[i,2]
    epcUnit <- paste0("(",epcUnit,")")
    npad2 <- 10 - nchar(epcUnit)
    if (grepl("\\*", epcUnit)){
      epcUnit <- gsub("\\(\\*", "\\*\\(", epcUnit)
      npad1 <- npad1 - 1
      npad2 <- npad2 + 1
    }
    epcDescription <- epc[i,1]
    if(epcUnit == "(flag)"){
      epcDescription <- strsplit(epcDescription, "; ")
      npad <- 22 - nchar(epcDescription[[1]][1])
      epcDescription <- paste0(
        epcDescription[[1]][1], 
        paste(rep(" ", npad), collapse = ""),
        epcDescription[[1]][2])
    }
    cat(epcValue, rep("", npad1))
    cat(epcUnit, rep("", npad2))
    cat(epcDescription)
  }
sink()
}