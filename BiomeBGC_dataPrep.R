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
  reqdPkgs = list("PredictiveEcology/SpaDES.core@box (>= 2.1.8.9013)", "ggplot2", "PredictiveEcology/BiomeBGCR@development"),
  parameters = bindrows(
    defineParameter("epcDataSource", "character", NA, NA, NA, 
                    paste("Three options:",
                          "1) The path(s) to the ecophysiological constants file(s).",
                          "There can be a single file, or one per study site.",
                          "2) One of the Biome-BGC default epc files (e.g., \"enf\")",
                          "A single default can be used or one per study site.",
                          "3) The source from which ecophysiological constants are extracted (TBD).")
                    ),
    defineParameter("metDataSource", "character", NA, NA, NA, 
                    paste("Two options:",
                          "1) The path(s) to the meteorological file(s).",
                          "There can be a single file, or one per study site.",
                          "2) Source of the meteorological data. PRISM? daymet? CHELSA (TBD)?")
                    ),
    defineParameter("CO2DataSource", "numeric|character", NA, NA, NA, 
                    paste("Three options:",
                          "1) A numeric to set a constant atmospheric CO2 concentraion (in ppm),",
                          "2) The path to a single file with annual variable CO2 concentration",
                          "3) The name of a data source to extract variable CO2 concentration (TBD).")
                    ),
    defineParameter("NDepositionLevel", "numeric", c(0, 2099, 0.0001), NA, NA, 
                    paste("A 3-number vector:",
                          "1) Keep nitrogen deposition level constant (0) or vary according to the time trajectory of CO2 mole fractions (1).",
                          "2) The reference year for N deposition (only used when N-deposition varies).",
                          "3) Industrial N deposition value.")
    ),
    defineParameter("dailyOutput", "numeric", c(1,2,3), NA, NA, 
                    paste("The indices of the daily output variable(s) requested.",
                          "There are >500 possible variables.")
    ),
    defineParameter("annualOutput", "numeric", c(1,2,3), NA, NA, 
                    paste("The indices of the daily output variable(s) requested.",
                          "There are >500 possible variables.")
    ),
    defineParameter("maxSpinupYears", "integer", 6000L, NA, NA, 
                    paste("The maximum number of simulation for a spinup run.")
    ),
    defineParameter("climateChangeOptions", "numberic", c(0, 0, 1, 1, 1), NA, NA, 
                    paste("Entries in the CLIMATE_CHANGE section of the ini file.",
                          "The entries are: offset for Tmax, offset for Tmin",
                          "multiplier for prcp, multiplier for vpd, and muliplier",
                          "for srad.")
    ),
    defineParameter("siteConstants", "numeric", c(1, 30, 50, 20, 977, 46.8, 0.2, 0.0001, 0.0008), NA, NA, 
                    paste("A vector with site information:",
                          "1: effective soil depth",
                          "2: sand percentage",
                          "3: silt percentage",
                          "4: clay percentage",
                          "5: site elevation in meters", 
                          "6: site latitude in decimal degrees",
                          "7: site shortwave albedo",
                          "8: annual rate of atmospheric nitrogen deposition",
                          "9: annual rate of symbiotic+asymbiotic nitrogen fixation")
    ),
    defineParameter("siteNames", "character", NA, NA, NA, 
                    paste("The names of the study sites. If not provided, a number from",
                          "1 to the number of study sites will be used.")
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
    expectsInput("studyArea", "SpatVector",
                desc = paste("Polygons to use as the study area. Must be supplied by the user.",
                             "One polygon per study site.")
                ),
    expectsInput("ecophysiologicalConstants", "data.frame", 
                 desc = paste("Ecophysiological constants. The first column is the",
                              "description of the constants, the second column is",
                              "the unit, and the following are the values of the constant.",
                              "There are either a single column if all study sites share",
                              "the same constant of one column per sites.")
                 ),
    expectsInput("meteorologicalData", "list", 
                 desc = paste("List of data.frames with the meteorological data",
                              "for each study site. The units are `deg C`",
                              "for Tmax, Tmin, and Tday, `cm` for prcp, `Pa` for",
                              "VPD, `W/m^2` for srad, and `s` for daylen.")
                 ),
    expectsInput("CO2concentration", "data.frame", 
                 desc = paste("CO2 concentration for each year.")
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

### template for your event1
prepareSpinupIni <- function(sim) {
  browser()
  # First read the ini template
  bbgcSpinup.ini <- iniRead(system.file("inputs/ini/template.ini", package = "BiomeBGCR"))
  
  # Set MET_INPUT section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "MET_INPUT", 1, 
                           file.path("inputs", "metdata", basename(basename(P(sim)$metDataSource)))
                           )
  # Set RESTART section
  # TODO: make sure that the 
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RESTART", c(1:4), c(0, 1, 0, 0))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RESTART", c(5,6), 
                           file.path("inputs", "restart",
                                    paste0(paste0(P(sim)$siteNames, ".restart"))
                           )
  )
  
  # Set TIME_DEFINE section
  nyear <- length(unique(sim$meteorologicalData$year))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 1, nyear) # number of year in the metdata
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 2, end(sim) - start(sim)) # number of simulation years
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 3, start(sim)) #first simulation year
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 4, 1) # 1 = spinup, 0 = normal simulation
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 5, P(sim)$maxSpinupYears) # max spinup years
  
  # Set CLIM_CHANGE section
  cc_values <- P(sim)$climateChangeOptions
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", c(1:5), cc_values)
  
  # Set CO2_CONTROL section
  # TODO: make sure that co2 is always constant for the spinup. If so, which year?
  # TODO: Get from external source
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "CO2_CONTROL", 1, 0) # Constant co2 concentration during spinup
  if (is.numeric(P(sim)$CO2DataSource)) {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "CO2_CONTROL", 2, 
                             P(sim)$CO2DataSource)
  } else {
    bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "CO2_CONTROL", 2, 
                             sim$CO2concentration[sim$CO2concentration$year == start(sim), "concentration"])
  }
  
  # Set SITE section
  # TODO: Get from an external source
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 1:9, P(sim)$siteConstants)
  
  # Set RAMP_NDEP section
  # TODO: Make sure that it is always constant during spinup
  # TODO: Get from external source
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RAMP_NDEP", 1:3, 
                           c(0, #0 = constant deposition
                             P(sim)$NDepositionLevel[2],
                             P(sim)$NDepositionLevel[3]
                           ))
  
  # Set EPC_FILE section
  if(P(sim)$epcDataSource %in% c("c3grass", "c4grass", "dbf", "dnf", "ebf", "enf", "shrub")){
    fileName <- paste0(P(sim)$epcDataSource, ".epc")
  } else {
    fileName <- basename(P(sim)$epcDataSource)
  }
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "EPC_FILE", 1, 
                           file.path("inputs", "epc", fileName))
  
  # Set W_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "W_STATE", 1:2, P(sim)$waterState)
  
  # Set C_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "C_STATE", 1:11, P(sim)$carbonState)
  
  # Set N_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "N_STATE", 1:2, P(sim)$nitrogenState)
  
  # Set OUTPUT_CONTROL section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "OUTPUT_CONTROL", 1:7, P(sim)$outputControl)
  
  # Set DAILY_OUTPUT section
  nDailyOuput <- length(P(sim)$dailyOutput)
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini,
                           "DAILY_OUTPUT",
                           1:(nDailyOuput + 1),
                           c(nDailyOuput, P(sim)$dailyOutput))
  
  # Set ANNUAL_OUTPUT section
  nAnnOuput <- length(P(sim)$annualOutput)
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini,
                           "ANNUAL_OUTPUT",
                           1:(nAnnOuput + 1),
                           c(nAnnOutput, P(sim)$annualOutput))
  
  return(invisible(sim))
}

### template for your event2
prepareIni <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere('studyArea', sim)) {
    stop("studyArea must be provided.")
  }
  
  
  if (!suppliedElsewhere('ecophysiologicalConstants', sim)) {
    if(is.na(P(sim)$epcDataSource)){
      stop("Either ecophysiologicalConstants or the parameter epcDataSource needs to be provided.")
    } else {
      dir.create(file.path(inputPath(sim), "epc"), showWarnings = FALSE)
      if(P(sim)$epcDataSource %in% c("c3grass", "c4grass", "dbf", "dnf", "ebf", "enf", "shrub")){
        fileName <- file.path(
          system.file("inputs", package = "BiomeBGCR"),
          "epc",
          paste0(P(sim)$epcDataSource, ".epc")
        )
      } else {
        fileName <- P(sim)$epcDataSource
      }
      newFileName <- file.path(inputPath(sim), "epc", basename(fileName))
      file.copy(fileName, 
                newFileName,
                overwrite = T
                )
      sim$ecophysiologicalConstants <- epcRead(newFileName)
    }
  } else if (!is.na(P(sim)$epcDataSource)) {
    message("Using provided ecophysiological constants, ignoring parameter epcDataSource.")
  }
  
  if (!suppliedElsewhere('meteorologicalData', sim)) {
    if(is.na(P(sim)$metDataSource)){
      stop("Either meteorologicalData or the parameter metDataSource needs to be provided.")
    } else {
      dir.create(file.path(inputPath(sim), "metdata"), showWarnings = FALSE)
      newFileName <- file.path(inputPath(sim), "metdata", basename(P(sim)$metDataSource))
      file.copy(P(sim)$metDataSource, 
                newFileName,
                overwrite = T
      )
      sim$meteorologicalData <- metRead(newFileName)
    }
  } else if (!is.na(P(sim)$metDataSource)) {
    message("Using provided meterological data, ignoring parameter metDataSource.")
  }
  
  if (!suppliedElsewhere('CO2concentration', sim)) {
    if(is.na(P(sim)$CO2DataSource)){
      stop("Either CO2concentration or the parameter CO2DataSource needs to be provided.")
    } else if (is.numeric(P(sim)$CO2DataSource)){
      sim$CO2concentration <- P(sim)$CO2DataSource
    } else {
      dir.create(file.path(inputPath(sim), "co2"), showWarnings = FALSE)
      newFileName <- file.path(inputPath(sim), "co2", basename(P(sim)$CO2DataSource))
      file.copy(P(sim)$CO2DataSource, 
                newFileName,
                overwrite = T
      )
      sim$CO2concentration <- CO2Read(newFileName)
    }
  } else if (!is.na(P(sim)$CO2DataSource)) {
    message("Using provided CO2 concentration data, ignoring parameter CO2DataSource.")
  }
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

