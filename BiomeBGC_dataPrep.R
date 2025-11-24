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
    defineParameter("dailyOutput", "numeric", c(20, 21, 38, 40, 42, 43, 44, 70, 509, 528, 620, 621, 622, 623, 624, 625, 626, 627, 636, 637, 638, 639, 579), NA, NA, 
                    paste("The indices of the daily output variable(s) requested.",
                          "There are >500 possible variables.")
    ),
    defineParameter("annualOutput", "numeric", c(545, 636, 637, 638, 639, 307), NA, NA, 
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
    defineParameter("waterState", "numeric", c(0, 0.5), NA, NA, 
                    paste("2-number vector for initial water conditions:",
                          "1: initial snowpack water content (kg/m2)",
                          "2: intial soil water content as a proportion of saturation (DIM)")
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
    defineParameter("nitrogenState", "numeric", c(0, 0), NA, NA, 
                    paste("2-number vector for initial nitrogen conditions:",
                          "1: litter nitrogen associated with labile litter carbon pool (kgN/m2)",
                          "2: soil mineral nitrogen pool (kgN/m2)")
    ),
    defineParameter("siteNames", "character", "site1", NA, NA, 
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

prepareSpinupIni <- function(sim) {
  # First read the ini template
  bbgcSpinup.ini <- iniRead(system.file("inputs/ini/template.ini", package = "BiomeBGCR"))
  
  # Set MET_INPUT section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "MET_INPUT", 1, 
                           file.path(inputPath(sim), "metdata", basename(basename(P(sim)$metDataSource)))
                           )
  # Set RESTART section
  # TODO: make sure that the 
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RESTART", c(1:4), c(0, 1, 0, 0))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini,
                           "RESTART",
                           c(5, 6),
                           file.path(inputPath(sim), "restart", paste0(P(sim)$siteNames, ".restart")))
  
  # Set TIME_DEFINE section
  nyear <- length(unique(sim$meteorologicalData$year))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 1, nyear) # number of year in the metdata
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 2, end(sim) - start(sim)) # number of simulation years
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 3, start(sim)) #first simulation year
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 4, 1) # 1 = spinup, 0 = normal simulation
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "TIME_DEFINE", 5, P(sim)$maxSpinupYears) # max spinup years
  
  # Set CLIM_CHANGE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "CLIM_CHANGE", c(1:5), P(sim)$climateChangeOptions)
  
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
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", c(1,5), P(sim)$siteConstants[c(1,5)])
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", c(2:4), P(sim)$siteConstants[c(2:4)])
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 6, P(sim)$siteConstants[6])
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", 7, P(sim)$siteConstants[7])
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "SITE", c(8:9), format(P(sim)$siteConstants[c(8:9)], scientific = FALSE, trim = TRUE))
  
  
  # Set RAMP_NDEP section
  # TODO: Make sure that it is always constant during spinup
  # TODO: Get from external source
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "RAMP_NDEP", 1:3, 
                           c(0, #0 = constant deposition
                             P(sim)$NDepositionLevel[2],
                             format(P(sim)$NDepositionLevel[3], scientific = FALSE, trim = TRUE)
                           ))
  
  # Set EPC_FILE section
  if(P(sim)$epcDataSource %in% c("c3grass", "c4grass", "dbf", "dnf", "ebf", "enf", "shrub")){
    fileName <- paste0(P(sim)$epcDataSource, ".epc")
  } else {
    fileName <- basename(P(sim)$epcDataSource)
  }
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "EPC_FILE", 1, 
                           file.path(inputPath(sim), "epc", fileName))
  
  # Set W_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "W_STATE", 1:2, P(sim)$waterState)
  
  # Set C_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "C_STATE", 1:11, P(sim)$carbonState)
  
  # Set N_STATE section
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "N_STATE", 1:2, P(sim)$nitrogenState)
  
  # Set OUTPUT_CONTROL section
  # TODO make sure this is what we want for the spinup
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "OUTPUT_CONTROL", 1, 
                           file.path(outputPath(sim), P(sim)$siteNames, "_spinup"))
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini, "OUTPUT_CONTROL", 2:6, 
                           c(0, # 1 = write daily output   0 = no daily output
                             0, # 1 = monthly avg of daily variables  0 = no monthly avg
                             0, # 1 = annual avg of daily variables   0 = no annual avg
                             1, # 1 = write annual output  0 = no annual output
                             1)) # for on-screen progress indicator
  
  # Set DAILY_OUTPUT section
  # TODO Get the comments right
  nDailyOuput <- length(P(sim)$dailyOutput)
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini,
                           "DAILY_OUTPUT",
                           1:(nDailyOuput + 1),
                           c(nDailyOuput, P(sim)$dailyOutput))
  bbgcSpinup.ini$DAILY_OUTPUT <- bbgcSpinup.ini$DAILY_OUTPUT[c(1:(2+nDailyOuput)),]
  
  # Set ANNUAL_OUTPUT section
  nAnnOuput <- length(P(sim)$annualOutput)
  bbgcSpinup.ini <- iniSet(bbgcSpinup.ini,
                           "ANNUAL_OUTPUT",
                           1:(nAnnOuput + 1),
                           c(nAnnOuput, P(sim)$annualOutput))
  bbgcSpinup.ini$ANNUAL_OUTPUT <- bbgcSpinup.ini$ANNUAL_OUTPUT[c(1:(2+nAnnOuput)),]
  
  # add to simList
  sim$bbgcSpinup.ini <- bbgcSpinup.ini
  
  return(invisible(sim))
}

prepareIni <- function(sim) {
  # Start from the spinup ini
  bbgc.ini <- sim$bbgcSpinup.ini
  
  # Change the RESTART section
  bbgc.ini <- iniSet(bbgc.ini, "RESTART", c(1:4), c(1, 0, 0, 0))
  bbgc.ini <- iniSet(bbgc.ini, "RESTART", 5, file.path(inputPath(sim), "restart", paste0(P(sim)$siteNames, ".restart")))
  bbgc.ini <- iniSet(bbgc.ini, "RESTART", 6, file.path(inputPath(sim), "restart", paste0(P(sim)$siteNames, "_out.restart")))
  
  # Change the TIME_DEFINE section
  nyear <- length(unique(sim$meteorologicalData$year))
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 1, nyear) # number of year in the metdata
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 2, end(sim) - start(sim)) # number of simulation years
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 3, start(sim)) #first simulation year
  bbgc.ini <- iniSet(bbgc.ini, "TIME_DEFINE", 4, 0) # 1 = spinup, 0 = normal simulation
  
  # Change the CO2 section
  if (!is.numeric(P(sim)$CO2DataSource)) {
    bbgc.ini <- iniSet(bbgc.ini, "CO2_CONTROL", 1, 1)
    bbgc.ini <- iniSet(bbgc.ini,
                       "CO2_CONTROL",
                       3,
                       file.path(inputPath(sim), "co2", basename(P(sim)$CO2DataSource)))
  }
  
  # Change the RAMP_NDEP section
  bbgc.ini <- iniSet(bbgc.ini, "RAMP_NDEP", 1, P(sim)$NDepositionLevel[1])
  
  # Change OUTPUT_CONTROL section
  bbgc.ini <- iniSet(bbgc.ini,
                     "OUTPUT_CONTROL",
                     1,
                     file.path(outputPath(sim), P(sim)$siteNames))
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

