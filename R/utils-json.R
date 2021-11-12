## Microbiome-specific json helpers

#' @import veupathUtils
writeAppResultsToJson <- function(appResults, pattern = NULL, dir = NULL, verbose = c(TRUE, FALSE)) {
  verbose <- matchArg(verbose)
  
  if (is.null(pattern)) pattern <- 'file'
  if (is.null(dir)) dir <- tempdir()
  
  # Simple formatting
  outJson <- getAppJson(appResults)
  
  outFileName <- basename(tempfile(pattern = pattern, tmpdir = dir, fileext = ".json"))
  write(outJson, outFileName)
  veupathUtils::logWithTime(paste('New json file written:', outFileName), verbose)
  
  return(outFileName)
}

#' @importFrom jsonlite unbox
#' @importFrom jsonlite toJSON
getAppJson <- function(appResults) {
  
  #### turn into some nice functions
  parameterSets <- lapply(appResults, function(dt) {return(attr(dt,'parameters'))})
  computations <- lapply(appResults, function(dt) {
    computation <- list()
    attr <- attributes(dt)
    
    ### Note - will need to ensure computedVariableDetails is written as an array of computedVariableDetails
    computation$computedVariableDetails <- attr$computedVariableDetails
    computation$computationDetails <- jsonlite::unbox(attr$computationDetails)
    
    # App-specific attributes
    computation$isCutoff <- jsonlite::unbox(attr$isCutoff)
    computation$pcoaVariance <- attr$pcoaVariance
    if ('recordVariable' %in% names(attr)) {
      computation$recordVariableDetails <- list('variableId' = veupathUtils::strSplit(attr$recordVariable,".", 4, 2),
                                                'entityId' = veupathUtils::strSplit(attr$recordVariable,".", 4, 1),
                                                'values' = dt[[attr$recordVariable]])
      dt[[attr$recordVariable]] <- NULL
    }
    
    # Set computation values
    computation$computedVariableDetails$values <- lapply(seq_along(dt), function(x) {return(dt[[x]])})
    
    return(computation)
  })
  
  outList <- list('computations' = computations,
                  'parameterSets' = parameterSets)
  
  outJson <- jsonlite::toJSON(outList)
  return(outJson)
}
