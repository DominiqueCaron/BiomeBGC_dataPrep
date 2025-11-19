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
                    paste("Three options:"
                          "1) The path(s) to the ecophysiological constants file(s).",
                          "There can be a single file, or one per study site.",
                          "2) One of the Biome-BGC default epc files (e.g., \"enf\")",
                          "A single default can be used or one per study site.",
                          "3) The source from which ecophysiological constants are extracted (TBD).")
                    ),
    defineParameter("metDataSource", "character", NA, NA, NA, 
                    paste("Two options:",
                          "1) The path(s) to the meteorological file(s).",
                          "There can be a single file, or one per study site."
                          "2) Source of the meteorological data. PRISM? daymet? CHELSA (TBD)?")
                    ),
    defineParameter("CO2DataSource", "numeric|character", NA, NA, NA, 
                    paste("Three options:",
                          "1) A numberic to set a constant atomospheric CO2 concentraion (in ppm),",
                          "2) The path to a single file with annual variable CO2 concentration",
                          "3) The name of a data source to extract variable CO2 concentration (TBD).")
                    ),
    defineParameter("NDepositionLevel", "numeric", 0, NA, NA, 
                    paste("Either `0` or a 2-number vector. If `0`, the nitrogen",
                          "deposition level remains constant.",
                          "If a 2-number vector, the nitrogen deposition levels",
                          "vary according to the time trajectory of CO2 mole fractions.",
                          "See Biome-BGC documentation.")
    ),
    defineParameter("dailyOutput", "numeric", c(1,2,3), NA, NA, 
                    paste("The indices of the daily output variable(s) requested.",
                          "There are >500 possible variables.")
    ),
    defineParameter("annualOutput", "numeric", c(1,2,3), NA, NA, 
                    paste("The indices of the daily output variable(s) requested.",
                          "There are >500 possible variables.")
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
    expectsInput("Co2concentration", "data.frame", 
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
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

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
    event1 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "BiomeBGC_dataPrep", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    event2 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "BiomeBGC_dataPrep", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #

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
  # do stuff for this event
  sampleData <- data.frame("TheSample" = sample(1:10, replace = TRUE))
  Plots(sampleData, fn = ggplotFn) # needs ggplot2

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
Event2 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

ggplotFn <- function(data, ...) {
  ggplot2::ggplot(data, ggplot2::aes(TheSample)) +
    ggplot2::geom_histogram(...)
}

