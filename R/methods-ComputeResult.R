#' Write Computed Variable Metadata
#'
#' This function, for any given ComputeResult, 
#' will return a filehandle where ComputedVariableMetadata
#' has been written in JSON format.
#' 
#' @param object ComputeResult 
#' @return filehandle where JSON representation of ComputedVariableMetadata can be found
#' @export
setGeneric("writeMeta",
  function(object, pattern = NULL, verbose = c(TRUE, FALSE)) standardGeneric("writeMeta"),
  signature = c("object")
)

#'@export 
setMethod("writeMeta", signature("ComputeResult"), function(object, pattern = NULL, verbose = c(TRUE, FALSE)) {
  verbose <- veupathUtils::matchArg(verbose)

  outJson <- veupathUtils::toJSON(object@computedVariableMetadata)

  if (is.null(pattern)) { 
    pattern <- object@name
    if (is.null(pattern)) {
      pattern <- 'file'
    } 
  }
  pattern <- paste0(pattern, '-meta-')

  outFileName <- basename(tempfile(pattern = pattern, tmpdir = tempdir(), fileext = ".json"))
  write(outJson, outFileName)
  veupathUtils::logWithTime(paste('New output file written:', outFileName), verbose)

  return(outFileName)
})

#' Write Compute Result Data
#'
#' This function, for any given ComputeResult, 
#' will return a filehandle where result data
#' has been written in tab delimited format.
#' @param object ComputeResult 
#' @return filehandle where tab delimited representation of result data can be found
#' @export
setGeneric("writeData",
  function(object, pattern = NULL, verbose = c(TRUE, FALSE)) standardGeneric("writeData"),
  signature = c("object")
)

#'@export 
setMethod("writeData", signature("ComputeResult"), function(object, pattern = NULL, verbose = c(TRUE, FALSE)) {
  verbose <- veupathUtils::matchArg(verbose)

  if (is.null(pattern)) { 
    pattern <- object@name
    if (is.null(pattern)) {
      pattern <- 'file'
    } 
  }
  pattern <- paste0(pattern, '-data-')

  outFileName <- basename(tempfile(pattern = pattern, tmpdir = tempdir(), fileext = ".tab"))
  data.table::fwrite(object@data, outFileName, sep = '\t', quote = FALSE)
  veupathUtils::logWithTime(paste('New output file written:', outFileName), verbose)

  return(outFileName)
})




#' Write Computed Variable Statistics
#'
#' This function, for any given ComputeResult, 
#' will return a filehandle where the statistics table
#' has been written in JSON format.
#' 
#' @param object ComputeResult with statistics
#' @return filehandle where JSON representation of ComputeResult statistics can be found
#' @export
setGeneric("writeStatistics",
  function(object, pattern = NULL, verbose = c(TRUE, FALSE)) standardGeneric("writeStatistics"),
  signature = c("object")
)

#'@export 
setMethod("writeStatistics", signature("ComputeResult"), function(object, pattern = NULL, verbose = c(TRUE, FALSE)) {
  verbose <- veupathUtils::matchArg(verbose)

  # Convert all to character but maintain structure 
  outObject <- data.frame(lapply(object@statistics, as.character))

  outJson <- jsonlite::toJSON(outObject)

  if (is.null(pattern)) { 
    pattern <- object@name
    if (is.null(pattern)) {
      pattern <- 'file'
    } 
  }
  pattern <- paste0(pattern, '-statistics-')

  outFileName <- basename(tempfile(pattern = pattern, tmpdir = tempdir(), fileext = ".json"))
  write(outJson, outFileName)
  veupathUtils::logWithTime(paste('New output file written:', outFileName), verbose)

  return(outFileName)
})