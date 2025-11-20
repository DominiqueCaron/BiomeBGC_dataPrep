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
  attr(epcDF, "decription") <- functionalType
  
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

}