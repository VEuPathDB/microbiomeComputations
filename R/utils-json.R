## Microbiome-specific json helpers

#' @import veupathUtils
writeAppResultsToJson <- function(appResults, recordIdColumn, pattern = NULL, dir = NULL, verbose = c(TRUE, FALSE)) {
  verbose <- matchArg(verbose)
  
  if (is.null(pattern)) pattern <- 'file'
  if (is.null(dir)) dir <- tempdir()
  
  # Simple formatting
  outJson <- getAppJson(appResults, recordIdColumn)
  
  outFileName <- basename(tempfile(pattern = pattern, tmpdir = dir, fileext = ".json"))
  write(outJson, outFileName)
  veupathUtils::logWithTime(paste('New json file written:', outFileName), verbose)
  
  return(outFileName)
}

#' @importFrom jsonlite unbox
#' @importFrom jsonlite toJSON
#' @import veupathUtils
getAppJson <- function(appResults, recordIdColumn) {
  
  #### turn into some nice functions
  parameterSets <- lapply(appResults, function(dt) {return(attr(dt,'parameters'))})
  recordVariableDetails <- list('variableId' = jsonlite::unbox(veupathUtils::strSplit(recordIdColumn,".", 4, 2)),
                                'entityId' = jsonlite::unbox(veupathUtils::strSplit(recordIdColumn,".", 4, 1)),
                                'values' = appResults[[1]][[recordIdColumn]])
  
  computations <- lapply(appResults, function(dt, recordIdColumn) {
    
    computation <- list()
    attr <- attributes(dt)
    
    # Set general computation-level attributes
    computation$computationDetails <- jsonlite::unbox(attr$computationDetails)
    
    # App/computation-specific attributes
    computation$isCutoff <- jsonlite::unbox(attr$isCutoff)
    computation$pcoaVariance <- attr$pcoaVariance

    # Set computation variable
    computation$computedVariable <- attr$computedVariable

    # Set computation values. Assuming only one computed variable for now.
    columnNames <- veupathUtils::toColNameOrNull(attr$computedVariable$computedVariableDetails)
    dt[[recordIdColumn]] <- NULL
    computation$computedVariable$computedVariableDetails$values <- lapply(seq_along(columnNames), function(x,dt) {return(dt[[x]])}, dt=dt)
    
    # Handle collection variable if present
    if (!is.null(attr$computedVariable$computedVariableMetadata$collectionVariable)) {
      computation$computedVariable$computedVariableMetadata$collectionVariable$collectionType <- jsonlite::unbox(as.character(attr$computedVariable$computedVariableMetadata$collectionVariable$collectionType))
    }

    # computation$computedVariable$values <- lapply(attr$computedVariable, function(computedVariable, dt) {
    #   return(computedVariable)
    # }, dt=dt)
    
    return(computation)
  }, recordIdColumn=recordIdColumn)
  
  outList <- list('computations' = computations,
                  'parameterSets' = parameterSets,
                  'recordVariableDetails' = recordVariableDetails)
  
  outJson <- jsonlite::toJSON(outList)
  return(outJson)
}
