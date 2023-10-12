#' Correlation
#'
#' This function returns correlation coefficients for variables in one dataset against variables in a second dataset
#' 
#' @param data1 first dataset. An AbundanceData object or data.table
#' @param data2 second dataset. A SampleMetadata object (if data1 is class AbundanceData) or a data.table
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
setGeneric("correlation",
  function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE), ...) standardGeneric("correlation"),
  signature = c("data1","data2")
)


#' Correlation
#'
#' This function returns correlation coefficients for all columns in one data table with all columns in a second data table.
#' 
#' @param data1 data.table with columns as variables. All columns must be numeric. One row per sample.
#' @param data2 data.table with columns as variables. All columns must be numeric. One row per sample. Will correlate all columns of data2 with all columns of data1.
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
setMethod("correlation", signature("data.table", "data.table"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {

  # Check that the number of rows match.
  if (!identical(nrow(data1), nrow(data2))) {
    stop("data1 and data2 must have the same number of rows.")
  }

  # Check that all values are numeric
  if (!identical(veupathUtils::findNumericCols(data1), names(data1))) { stop("All columns in data1 must be numeric.")}
  if (!identical(veupathUtils::findNumericCols(data2), names(data2))) { stop("All columns in data2 must be numeric.")}


  ## Compute correlation
  # Resulting data table has column "rn" = row names of the correlation matrix (so taxa names), and the other
  # column names are the vars from sample metadata that we use
  corrResult <- data.table::as.data.table(cor(data1, data2, method = method), keep.rownames = T)

  veupathUtils::logWithTime(paste0('Completed correlation with method=', method,'. Formatting results.'), verbose)


  ## Format results
  meltedCorrResult <- melt(corrResult, id.vars=c('rn'))
  formattedCorrResult <- data.frame(
    data1 = meltedCorrResult[['rn']],
    data2 = meltedCorrResult[['variable']],
    correlationCoef = meltedCorrResult[['value']]
  )

  return(formattedCorrResult)

})

#' Correlation of abundance data and metadata
#' 
#' This function returns the correlation of all columns of an AbundanceData object with appropriate columns of a SampleMetadata object.
#' 
#' @param data1 AbundanceData object. Will correlate abundance variables with specified variables in data2
#' @param data2 SampleMetadata object. Will assume that all non-ID columns of metadata are appropriate for correlation.
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @export
setMethod("correlation", signature("AbundanceData", "SampleMetadata"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {

  # Many of our checks betwee smaple metadata and data live in the object check function of the abudnance class. 
  # So we can use those (runnning validObject) to do a bunch of checks, then access the data.
  if (!is.null(data1@sampleMetadata)) {
    warning("Replacing sampleMetadata in data1 with data2 sampleMetadata for validation")
  }
  data1@sampleMetadata <- data2
  validObject(data1) # Checks that record IDs match


  df1 <- getAbundances(data1, FALSE) # use the imputeZero flag if it exists
  recordIdColumn <- data1@recordIdColumn
  ancestorIdColumns <- data1@ancestorIdColumns
  allIdColumns <- c(recordIdColumn, ancestorIdColumns)
  df2 <- getSampleMetadata(data1) # returns a data.table



  ## Initialize and check inputs
  method <- veupathUtils::matchArg(method)
  verbose <- veupathUtils::matchArg(verbose)


  computeMessage <- ''
  veupathUtils::logWithTime(paste("Received df table with", nrow(df1), "samples and", (ncol(df1)-1), "taxa."), verbose)


  ## Compute correlation
  # Result has columns "data1", "data2", and correlationCoef
  corrResult <- correlation(df1[, -..allIdColumns], df2[, -..allIdColumns], method = method, verbose = verbose)


  ## Format results
  # Construct the ComputeResult
  result <- new("CorrelationComputeResult")
  result@name <- 'correlation'
  result@recordIdColumn <- recordIdColumn
  result@ancestorIdColumns <- ancestorIdColumns
  result@statistics <- corrResult
  result@parameters <- paste0('method = ', method)
  result@data1Metadata <- "assay"
  result@data2Metadata <- "sampleMetadata"


  # The resulting data should contain only the samples actually used.
  result@data <- df1[, ..allIdColumns]
  names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))


  validObject(result)
  veupathUtils::logWithTime(paste('Correlation computation completed with parameters recordIdColumn=', recordIdColumn, ', method = ', method), verbose)
  
  return(result)
})