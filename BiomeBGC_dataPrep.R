## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "BiomeBGC_dataPrep",
  description = "Prepares inputs to run Biome-BGC.",
  keywords = "",
  authors = c(
    person("Dominique", "Caron", email = "dominique.caron@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = "aut"),
    person("Céline", "Boisvenue", email = "celine.boisvenue@nrcan-rncan.gc.ca", role = "ctb")
  ),
  childModules = character(0),
  version = list(BiomeBGC_dataPrep = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "BiomeBGC_dataPrep.Rmd"),
  reqdPkgs = list("PredictiveEcology/SpaDES.core@box (>= 2.1.8.9013)", "ggplot2", "PredictiveEcology/LandR@development",
                  "PredictiveEcology/BiomeBGCR@development", "elevatr", "terra", "rvest", "data.table",
                  "RNCan/BioSimClient_R", "geosphere"),
  parameters = bindrows(
    defineParameter("annualOutput", "numeric", c(545, 636, 637, 638, 639, 307), NA, NA, 
                    paste("The indices of the daily output variable(s) requested.",
                          "There are >500 possible variables.")
    ),
    defineParameter("carbonState", "numeric", c(0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), NA, NA, 
                    paste("11-number vector for initial carbon conditions:",
                          "1: peak leaf carbon to be attained during the first simulation year (kgC/m2)",
                          "2: peak stem carbon to be attained during the first year (kgC/m2)",
                          "3: initial coarse woody debris carbon (dead trees, standing or fallen) (kgC/m2)",
                          "4: initial litter carbon, labile pool (kgC/m2)",
                          "5: initial litter carbon, unshielded cellulose pool (kgC/m2)",
                          "6: initial litter carbon, shielded cellulose pool (kgC/m2)",
                          "7: initial litter carbon, lignin pool (kgC/m2)",
                          "8: soil carbon, fast pool (kgC/m2)",
                          "9: soil carbon, medium pool (kgC/m2)",
                          "10: soil carbon, slow pool (kgC/m2)",
                          "11: soil carbon, slowest pool (kgC/m2)"
                    )
    ),
    defineParameter("climModel", "character", "RCM4", NA, NA, 
                    paste("A climatic model to extract meteorological data.",
                          "Either 'RCM4', 'GCM4', or 'Hadley'.")
    ),
    defineParameter("climateChangeOptions", "numeric", c(0, 0, 1, 1, 1), NA, NA, 
                    paste("Entries in the CLIMATE_CHANGE section of the ini file.",
                          "The entries are: offset for Tmax, offset for Tmin",
                          "multiplier for prcp, multiplier for vpd, and muliplier",
                          "for srad.")
    ),
    defineParameter("co2scenario", "character", "RCP45", NA, NA, 
                    paste("An representative concentration pathway for the co2",
                          "concentration trajectories and meteorological data.",
                          "Either 'RCP45' or 'RCP85'.")
    ),
    defineParameter("dailyOutput", "numeric", c(20, 21, 38, 40, 42, 43, 44, 70, 509, 528, 620, 621, 622, 623, 624, 625, 626, 627, 636, 637, 638, 639, 579), NA, NA, 
                    paste("The indices of the daily output variable(s) requested.",
                          "There are >500 possible variables.")
    ),
    defineParameter("maxSpinupYears", "integer", 6000L, NA, NA, 
                    paste("The maximum number of simulation for a spinup run.")
    ),
    defineParameter("metSpinupYears", "numeric", 40, NA, NA, 
                    paste("The number of years used for the spinup.")
    ),
    defineParameter("NDepositionLevel", "numeric", c(1, NA, NA), NA, NA, 
                    paste("A 3-number vector:",
                          "1) Keep nitrogen deposition level constant (0) or vary according to the time trajectory of CO2 mole fractions (1).",
                          "2) The reference year for N deposition (only used when N-deposition varies).",
                          "3) Industrial N deposition value.")
    ),
    defineParameter("nitrogenState", "numeric", c(0, 0), NA, NA, 
                    paste("2-number vector for initial nitrogen conditions:",
                          "1: litter nitrogen associated with labile litter carbon pool (kgN/m2)",
                          "2: soil mineral nitrogen pool (kgN/m2)")
    ),
    defineParameter("siteConstants", "numeric", c(1, NA, NA, NA, NA, NA, NA, NA, NA), NA, NA, 
                    paste("A vector with site information:",
                          "1: effective soil depth",
                          "2: sand percentage",
                          "3: silt percentage",
                          "4: clay percentage",
                          "5: site elevation in meters", 
                          "6: site latitude in decimal degrees",
                          "7: site shortwave albedo",
                          "8: annual rate of atmospheric nitrogen deposition",
                          "9: annual rate of symbiotic+asymbiotic nitrogen fixation",
                          "The non-na constants will be retrieved in various sources.")
    ),
    defineParameter("siteNames", "character", "site1", NA, NA, 
                    paste("The names of the study sites.")
    ),
    defineParameter("waterState", "numeric", c(NA, 0.5), NA, NA, 
                    paste("2-number vector for initial water conditions:",
                          "1: initial snowpack water content (kg/m2)",
                          "2: intial soil water content as a proportion of saturation (DIM)",
                          "If set to NA, the initial snowpack water content will be retrieved by an external source.")
    ),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                    "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("climatePolygons", "SpatVector", 
                 desc = paste("Polygons of homogeneous climate. By default,",
                              "ecodistricts are used.")
    ),
    expectsInput("CO2concentration", "data.frame", 
                 desc = paste("CO2 concentration for each year.")
    ),
    expectsInput("dominantSpecies", "SpatRaster",
                 desc = paste(
                   "A raster with the leading tree species (speciesId).",
                   "Use to determine the ecophysiological constants.",
                   "By default, the NTEMS dominant tree species layer for the starting year is used.")
    ), 
    expectsInput("ecophysiologicalConstants", "data.frame", 
                 desc = paste(
                   "Ecophysiological constants. Columns are species, genus, and",
                   "the ecological constants (see template). By default,",
                   "ecophysiologicalConstants are built from White et al., 2000.")
    ),
    expectsInput("elevation", "SpatRaster",
                 desc = paste(
                   "An elevation (m) raster. By default, the raster will be extracted from",
                   "AWS Terrain Tiles with the elevatr package."
                 ),
                 sourceURL = "https://registry.opendata.aws/terrain-tiles/"
    ),
    expectsInput("meteorologicalData", "list", 
                 desc = paste("List of data.frames with the meteorological data",
                              "for each climate polygons. The units are `deg C`",
                              "for Tmax, Tmin, and Tday, `cm` for prcp, `Pa` for",
                              "VPD, `W/m^2` for srad, and `s` for daylen.")
    ),
    expectsInput("Ndeposition", "SpatRaster",
                 desc = paste(
                   "Raster(s) of total atmospheric N deposition (kgN/m2/yr).",
                   "If N deposition is variable two layers need to be provided,", 
                   "one for N deposition at the start, of the simulation and a",
                   "second for N deposition at another timestep. The layer name", 
                   "of the second raster needs to be the year of the data."
                 ),
                 sourceURL = "https://www.nature.com/articles/s41467-024-55606-y"
    ), 
    expectsInput("NFixationRates", "SpatRaster",
                 desc = paste(
                   "Raster of annual rate of symbiotic + asymbiotic nitrogen fixation (kgN/m2/yr)."
                 ),
                 sourceURL = "https://www.sciencebase.gov/catalog/item/66a97480d34e07a119db3a37"
    ), 
    expectsInput("pixelGroupMap", "SpatRaster", 
                 desc = paste("")
    ),
    expectsInput("pixelGroupParameters", "data.frame", 
                 desc = paste("")
    ),
    expectsInput("rasterToMatch", "SpatRaster", 
                 desc = paste("A raster defining the extent, resolution, projection of the",
                              "study area.")
    ),
    expectsInput("snowpackWaterContent", "SpatRaster",
                 desc = paste(
                   "Initial snowpack water content (kg/m2)."
                 ),
                 sourceURL = "https://climate-scenarios.canada.ca/?page=blended-snow-data"
    ), 
    expectsInput("soilTextures", "SpatRaster",
                 desc = paste(
                   "A raster stack with layers representing the % of 'Sand', 'Silt', and 'Clay'.",
                   "The across-layers sum needs to equal to 1 for each pixels."
                 ),
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/index.html"
    ), 
    expectsInput("sppEquiv", "data.frame",
                 desc = paste(
                   "A data frame to link the leading species map to ecophysiological constants.",
                   "The columns are speciesId, species, genus, functional plant type.")
    ), 
    expectsInput("studyArea", "SpatVector",
                 desc = paste("Polygons to use as the study area. Must be supplied by the user.",
                              "One polygon per study site.")
    )
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "bbgcSpinup.ini", objectClass = "character", desc = paste("Biome-BGC initialization files for the spinup.", 
                                                                                         "Path to the .ini files (one path per site/scenario).")),
    createsOutput(objectName = "bbgc.ini", objectClass = "character", desc =  paste("Biome-BGC initialization files.", 
                                                                                    "Path to the .ini files (one path per site/scenario)."))
  )
))

doEvent.BiomeBGC_dataPrep = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      
      sim <- prepareSpinupIni(sim)
      
      sim <- prepareIni(sim)
      
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "BiomeBGC_dataPrep", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "BiomeBGC_dataPrep", "save")
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      
      plotFun(sim) # example of a plotting function
      # schedule future event(s)
      
      # e.g.,
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "BiomeBGC_dataPrep", "plot")
      
      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function
      
      # schedule future event(s)
      
      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "BiomeBGC_dataPrep", "save")
      
      # ! ----- STOP EDITING ----- ! #
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}


### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

prepareSpinupIni <- function(sim) {
  # First read the ini template
  bbgcSpinup.ini <- iniRead(system.file("inputs/ini/template.ini", package = "BiomeBGCR"))
  
  ## Set MET_INPUT section
  fileName <- tolower(paste0(
    P(sim)$siteName,
    "_",
    P(sim)$climModel,
    P(sim)$co2scenario,
    "_spinup.mtc43"
  ))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "MET_INPUT", 1, 
                           file.path("inputs", "metdata", fileName))
  ## Set RESTART section
  # TODO: make sure that the 
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RESTART", c(1:4), c(0, 1, 0, 0))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini,
                           "RESTART",
                           c(5, 6),
                           file.path("inputs", "restart", paste0(P(sim)$siteNames, ".restart")))
  
  ## Set TIME_DEFINE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 1, P(sim)$metSpinupYears) # number of year in the metdata
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 2, P(sim)$metSpinupYears) # number of simulation years
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 3, start(sim)) #first simulation year
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 4, 1) # 1 = spinup, 0 = normal simulation
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 5, P(sim)$maxSpinupYears) # max spinup years
  
  ## Set CLIM_CHANGE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "CLIM_CHANGE", c(1:5), P(sim)$climateChangeOptions)
  
  ## Set CO2_CONTROL section
  # TODO: make sure that co2 is always constant for the spinup. If so, which year?
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "CO2_CONTROL", 1, 0) # Constant co2 concentration during spinup
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "CO2_CONTROL", 2, sim$CO2concentration[sim$CO2concentration$year == start(sim), "co2_ppm"])
  
  ## Set SITE section
  # For each check if NA, if it is get the value from different sources
  # Rooting zone soil depth, if NA set to 1m
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 1, 
                           ifelse(is.na(P(sim)$siteConstants[1]), 
                                  1, 
                                  P(sim)$siteConstants[1]))
  
  # Soil texture: % of sand, % of silt, % of clay
  if(is.na(P(sim)$siteConstants[2])){
    soilTexture <- extract(sim$soilTexture, sim$studyArea) |> round(digits = 1)
    # Make sure that sum == 100 after rounding.
    soilTexture[, "silt"] <- 100 - (soilTexture[, "sand"] + soilTexture[, "clay"])
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", c(2:4),
                             soilTexture[,-1])
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", c(2:4),
                             P(sim)$siteConstants[c(2:4)])
  }
  
  # Elevation
  if(is.na(P(sim)$siteConstants[5])){
    elevation <- extract(sim$elevation, sim$studyArea)[,2]
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 5, elevation)
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 5, P(sim)$siteConstants[5])
  }
  
  # Latitude
  if (is.na(P(sim)$siteConstants[6])) {
    latitude <- round(crds(project(sim$studyArea, "+proj=longlat +ellps=WGS84 +datum=WGS84"))[2], 2)
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 6, latitude)
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 6, P(sim)$siteConstants[6])
  }
  
  # Site shortwave albedo
  if (is.na(P(sim)$siteConstants[7])) {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 7, sim$shortwaveAlbedo)
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 7, P(sim)$siteConstants[7])
  }
  
  # wet+dry atmospheric deposition of N
  if (is.na(P(sim)$siteConstants[8])) {
    Ndeposition <- extract(sim$Ndeposition[[1]], sim$studyArea)[,2] |> round(digits = 5)
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 8, format(Ndeposition, scientific = FALSE, trim = TRUE))
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 8, P(sim)$siteConstants[8])
  }
  
  # symbiotic+asymbiotic fixation of N
  if (is.na(P(sim)$siteConstants[9])) {
    NfixationRate <- extract(sim$NfixationRates, sim$studyArea)[,2] |> round(digits = 4)
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini,
                             "SITE",
                             9,
                             format(NfixationRate, scientific = FALSE, trim = TRUE))
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 9, P(sim)$siteConstants[9])
  }
  
  # Set RAMP_NDEP section
  # TODO: Make sure that it is always constant during spinup
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RAMP_NDEP", 1, 0)
  if(P(sim)$NDeposition[1] == 1 & is.na(P(sim)$NDeposition[2])){
    year2 <- names(sim$Ndeposition[[2]])
    Ndeposition2 <- extract(sim$Ndeposition[[2]], sim$studyArea)[,2] |> round(digits = 5)
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RAMP_NDEP", 2, year2)
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RAMP_NDEP", 3, format(Ndeposition2, scientific = FALSE, trim = TRUE))
  } else if (P(sim)$NDeposition[1] == 1 & !is.na(P(sim)$NDeposition[2])){
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RAMP_NDEP", c(2,3), P(sim)$NDeposition[c(2,3)])
  }
  
  # Set EPC_FILE section
  # extract the correct dominant species
  dominantSpecies <- extract(sim$dominantSpecies, sim$studyArea)[-1]
  dominantSpecies <- unique(names(dominantSpecies)[dominantSpecies==1])
  if(length(dominantSpecies) == 0){
    stop("Study area is not a forested treed-area.")
  }
  dominantSpecies <- sim$sppEquiv$species[sim$sppEquiv$speciesId == dominantSpecies]
  # set filename
  fileName <- tolower(paste0(gsub(" ", "", dominantSpecies), ".epc"))
  # set in ini file
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "EPC_FILE", 1, 
                           file.path("inputs", "epc", fileName))
  
  # Set W_STATE section
  if(is.na(P(sim)$waterState[1])){
    snowpackWaterContent <- extract(sim$snowpackWaterContent, sim$studyArea)[, 2] |> round(digits = 1)
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "W_STATE", 1, snowpackWaterContent)
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "W_STATE", 1, P(sim)$waterState[1])
  }
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "W_STATE", 2, P(sim)$waterState[2])
  
  # Set C_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "C_STATE", 1:11, P(sim)$carbonState)
  
  # Set N_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "N_STATE", 1:2, P(sim)$nitrogenState)
  
  # Set OUTPUT_CONTROL section
  # TODO make sure this is what we want for the spinup
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "OUTPUT_CONTROL", 1, 
                           file.path("outputs", paste0(P(sim)$siteNames, "_spinup")))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "OUTPUT_CONTROL", 2:6, 
                           c(0, # 1 = write daily output   0 = no daily output
                             0, # 1 = monthly avg of daily variables  0 = no monthly avg
                             0, # 1 = annual avg of daily variables   0 = no annual avg
                             1, # 1 = write annual output  0 = no annual output
                             1)) # for on-screen progress indicator
  
  # Set DAILY_OUTPUT section
  nDailyOutput <- length(P(sim)$dailyOutput)
  # Rewrite, it is easier given that the size of the section varies
  bbgcSpinup.ini[["DAILY_OUTPUT"]] <- rbind(bbgcSpinup.ini[["DAILY_OUTPUT"]][1, ],
                                            c(nDailyOutput, "(int)", "number of daily variables to output"))
  bbgcSpinup.ini[["DAILY_OUTPUT"]] <- rbind(bbgcSpinup.ini[["DAILY_OUTPUT"]],
                                            data.frame(
                                              value = P(sim)$dailyOutput,
                                              unit = c(0:(nDailyOutput-1)),
                                              comment = getOutputDescription(P(sim)$dailyOutput)
                                            ))
  
  # Set ANNUAL_OUTPUT section
  nAnnOutput <- length(P(sim)$annualOutput)
  # Rewrite, it is easier given that the size of the section varies
  bbgcSpinup.ini[["ANNUAL_OUTPUT"]] <- rbind(bbgcSpinup.ini[["ANNUAL_OUTPUT"]][1, ],
                                             c(nAnnOutput, "(int)", "number of annual output variables"))
  bbgcSpinup.ini[["ANNUAL_OUTPUT"]] <- rbind(bbgcSpinup.ini[["ANNUAL_OUTPUT"]],
                                             data.frame(
                                               value = P(sim)$annualOutput,
                                               unit = c(0:(nAnnOutput-1)),
                                               comment = getOutputDescription(P(sim)$annualOutput)
                                             ))
  
  # add to simList
  sim$bbgcSpinup.ini <- bbgcSpinup.ini
  
  return(invisible(sim))
}

prepareIni <- function(sim) {
  # Start from the spinup ini
  bbgc.ini <- sim$bbgcSpinup.ini
  
  ## Set MET_INPUT section
  fileName <- tolower(paste0(
    P(sim)$siteName,
    "_",
    P(sim)$climModel,
    P(sim)$co2scenario,
    "_",
    start(sim),
    end(sim),
    ".mtc43"
  ))
  bbgc.ini <- iniSet(bbgc.ini,
                     "MET_INPUT",
                     1,
                     file.path("inputs", "metdata", fileName))
  
  # Change the RESTART section
  bbgc.ini <- iniSet(bbgc.ini, "RESTART", c(1:4), c(1, 0, 0, 0))
  bbgc.ini <- iniSet(bbgc.ini, "RESTART", 5, file.path("inputs", "restart", paste0(P(sim)$siteNames, ".restart")))
  bbgc.ini <- iniSet(bbgc.ini, "RESTART", 6, file.path("inputs", "restart", paste0(P(sim)$siteNames, "_out.restart")))
  
  # Change the TIME_DEFINE section
  nyear <- length(unique(sim$meteorologicalData$year))
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 1, nyear) # number of year in the metdata
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 2, end(sim) - start(sim) + 1) # number of simulation years
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 3, start(sim)) #first simulation year
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 4, 0) # 1 = spinup, 0 = normal simulation
  
  # Change the CO2 section
  fileName <- paste("co2", start(sim), end(sim), paste0(P(sim)$co2scenario, ".txt"), sep = "_")
  bbgc.ini <- iniSet(bbgc.ini, "CO2_CONTROL", 1, 1)
  bbgc.ini <- iniSet(bbgc.ini,
                     "CO2_CONTROL",
                     3,
                     file.path("inputs", "co2", fileName))
  
  # Change the RAMP_NDEP section
  bbgc.ini <- iniSet(bbgc.ini, "RAMP_NDEP", 1, P(sim)$NDepositionLevel[1])
  
  # Change OUTPUT_CONTROL section
  bbgc.ini <- iniSet(bbgc.ini,
                     "OUTPUT_CONTROL",
                     1,
                     file.path("outputs", P(sim)$siteNames))
  bbgc.ini <- iniSet(bbgc.ini, "OUTPUT_CONTROL", 2:6, c(
    1, # 1 = write daily output   0 = no daily output
    1, # 1 = monthly avg of daily variables  0 = no monthly avg
    1, # 1 = annual avg of daily variables   0 = no annual avg
    1, # 1 = write annual output  0 = no annual output
    1  # for on-screen progress indicator
  )) 
  
  sim$bbgc.ini <- bbgc.ini
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  # Study area needs to be either points or a polygon
  if (!suppliedElsewhere('studyArea', sim)) {
    stop("studyArea must be provided.")
  }
  
  # Climate polygons: Climate is assumed to be homogeneous within polygons
  # Default source: Canadian ecodistrict
  if (!suppliedElsewhere('climatePolygons', sim)) {
    
    sim$climatePolygons <- prepInputs(
      targetFile = "ecodistricts.shp",
      url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
      destinationPath = dPath,
      to = sim$rasterToMatch,
      fun = "terra::vect"
    ) |> Cache()
    
  }
  
  # Dominant species layer
  # Default source: NTEMS dominant species layer for the 1st year of simulation
  if (!suppliedElsewhere('dominantSpecies', sim)) {
    yearToUse <- max(min(start(sim), 2022), 1984)
    if (yearToUse != start(sim)){
      message("NTEMS dominant species layer is not available for ", start(sim), 
              ", using layer for ", yearToUse)
    }
    sim$dominantSpecies <- prepNTEMSDominantSpecies(
      year = yearToUse,
      destinationPath = dPath,
      cropTo = sim$rasterToMatch,
      projectTo = sim$rasterToMatch,
      maskTo = sim$studyArea
    ) |> Cache()
  }
  
  # Table to link the dominant species to traits of White et al., 2000
  if (!suppliedElsewhere('sppEquiv', sim)) {
    sppEquiv <- LandR::sppEquivalencies_CA[, c("LandR", "Latin_full", "Type")]
    sppEquiv <- sppEquiv[sppEquiv$Latin_full != "",]
    sppEquiv$genus <- sub(" .*", "", sppEquiv$Latin_full)
    sppEquiv$PFT <- ifelse(sppEquiv$Type == "Conifer", "enf", "dbf")
    sppEquiv <- data.frame(
      speciesId = sppEquiv$LandR,
      species = sppEquiv$Latin_full,
      genus = sppEquiv$genus,
      PFT = sppEquiv$PFT
    )
    speciesIdInStudyArea <- unique(values(sim$dominantSpecies)) |> na.omit()
    speciesInStudyArea <- cats(sim$dominantSpecies)[[1]][speciesIdInStudyArea,"category"]
    sim$sppEquiv <- sppEquiv[match(speciesInStudyArea, sppEquiv$speciesId),]
  }
  
  # Table of ecophysiological constants of each species Id
  # Default source: White et al., 2000 https://doi.org/10.1175/1087-3562(2000)004<0003:PASAOT>2.0.CO;2
  if (!suppliedElsewhere('ecophysiologicalConstants', sim)) {
    
    sim$ecophysiologicalConstants <- prepWhite2010EPC(
      url = "https://drive.google.com/file/d/1xVNwNenJRXtBKDTmpxMutS7Ipd5NSeWK/view?usp=drive_link",
      sppEquiv = sim$sppEquiv,
      destinationPath = dPath
    ) |> Cache()
    
  }
  
  # Soil texture (% sand, % silt, % clay)
  # Default sources: CanSIS Soil Landscape Grids of Canada
  # Using the 15-30cm depth layer
  if (!suppliedElsewhere('soilTexture', sim)) {
    sim$soilTexture <- prepSoilTexture(
      destinationPath = dPath,
      to = sim$rasterToMatch
    ) |> Cache()
  }
  
  # Elevation raster
  # Default source: Amazon Web Services Terrain Tiles
  if (!suppliedElsewhere('elevation', sim)) {
    elevation <-  get_elev_raster(
      locations = sf::st_as_sf(sim$studyArea),
      z = 10
    ) |> Cache()
    elevation <- postProcessTo(
      rast(elevation),
      to = sim$rasterToMatch
    ) |> Cache()
    sim$elevation <- signif(elevation, 3)
  }
  
  # Total N deposition
  # Default source: Zhu et al., 2025: https://doi.org/10.1038/s41467-024-55606-y
  if (!suppliedElsewhere('Ndeposition', sim)) {
    year1 <- max(start(sim), 2008)
    year2 <- min(end(sim), 2020)
    sim$Ndeposition <-  prepNdeposition(
      destinationPath = dPath,
      to = sim$rasterToMatch,
      year1 = year1,
      year2 = year2
    ) |> Cache()
  }
  
  # Total N fixation rates
  # Default source Reis Ely et al., 2025: https://doi.org/10.1038/s41597-025-05131-4
  if (!suppliedElsewhere('NFixationRates', sim)) {
    sim$NfixationRates <- prepInputs(
      targetFile = "BNF_total_central_1.tif",
      overwrite = TRUE,
      url = "https://www.sciencebase.gov/catalog/file/get/66b53cc6d34eebcf8bb3850e?f=__disk__67%2Fdf%2F6a%2F67df6a59f896d547205ddb20da99ec72db7a6b10",
      destinationPath = dPath,
      to = sim$rasterToMatch,
      fun = "terra::rast"
    ) |> Cache()
    sim$NfixationRates <- round(sim$NfixationRates, digits = 1)/10000 # convert from kg/ha/yr to kg/m2/yr
  }

  # Initial snowpack water content
  # Default source: ECC snow water equivalent (SWE) over the Northern Hemisphere
  # Methods: Mudryk et al., 2015: https://doi.org/10.1175/JCLI-D-15-0229.1
  if (!suppliedElsewhere('snowpackWaterContent', sim)) {
    
    # data is available for 1981-2020
    yearToUse <- min(max(start(sim), 1981), 2020)
    if(yearToUse != start(sim)){
      message("Snowpack water content data is not available for ",
              start(sim), ", using data of ", yearToUse)
    }
    
    # Get data
    sim$snowpackWaterContent <- prepInputs(
      targetFile = "swe_monthly_mm_1981-2020.nc",
      url = "https://climate-scenarios.canada.ca/files/blended_snow_2024/swe_monthly_mm_1981-2020.zip",
      fun = "terra::rast",
      destinationPath = dPath,
      to = sim$rasterToMatch,
    ) |> Cache()
    
    # We use the average for January of the first year.
    layerToKeep <- which(terra::time(sim$snowpackWaterContent) == paste(yearToUse, "01", "16", sep = "-"))
    sim$snowpackWaterContent <- round(sim$snowpackWaterContent[[layerToKeep]], 0)
    terra::units(sim$snowpackWaterContent) <- "kg/m^2"
  }
  
  # Shortwave Albedo
  # Default source: Based on SCANFI landcover and albedo of land cover type in 
  # Gao et al., 2005 (https://doi.org/10.1029/2004JD005190)
  if (!suppliedElsewhere('shortwaveAlbedo', sim)) {
    lcc <- prepInputs(
      url = "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/SCANFI/v1/SCANFI_att_nfiLandCover_SW_2020_v1.2.tif",
      destinationPath = dPath,
      to = sim$rasterToMatch,
      method = "near"
    ) |> Cache()
    albedoTable <- rvestAlbedoTable()
    sim$shortwaveAlbedo <- lccToAlbedo(lcc, albedoTable, sim$rasterToMatch)
  }
  
  # Meteorological data
  # Default source: BioSIM
  if (!suppliedElsewhere('meteorologicalData', sim)) {
    
    sim$meteorologicalData <- prepClimate(
      climatePolygons = sim$climatePolygons,
      siteName = P(sim)$siteNames,
      firstYear = start(sim),
      lastYear = end(sim),
      lastSpinupYear = start(sim) + P(sim)$metSpinupYears,
      scenario = P(sim)$co2scenario,
      climModel = P(sim)$climModel,
      destinationPath= dPath
    )|> Cache()
  }
  
  # CO2 atmospheric concentration
  # Default source: Emission data for the different RCPs downloaded from
  # Folini et al., 2025: https://doi.org/10.1093/restud/rdae011
  if (!suppliedElsewhere('CO2concentration', sim)) {
    
    sim$CO2concentration <- prepCo2Concentration(
      firstYear = start(sim),
      lastYear = end(sim),
      scenario = P(sim)$co2scenario,
      destinationPath= dPath
    ) |> Cache()
    
  }
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

