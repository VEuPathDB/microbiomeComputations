#' Correlation
#'
#' This function returns correlation coefficients for variables in one dataset against variables in a second dataset
#' 
#' @param data1 first dataset. An AbundanceData object or data.table
#' @param data2 second dataset. A SampleMetadata object (if data1 is class AbundanceData) or a data.table
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return data.frame with correlation coefficients or a ComputeResult object
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
#' @return data.frame with correlation coefficients
setMethod("correlation", signature("data.table", "data.table"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {

  # Check that the number of rows match.
  if (!identical(nrow(data1), nrow(data2))) {
    stop("data1 and data2 must have the same number of rows.")
  }

  # Check that all values are numeric
  if (!identical(veupathUtils::findNumericCols(data1), names(data1))) { stop("All columns in data1 must be numeric.")}
  if (!identical(veupathUtils::findNumericCols(data2), names(data2))) { stop("All columns in data2 must be numeric.")}


  ## Compute correlation
  # Resulting data table has column "rn" = row names of the correlation matrix (so taxa names, for example), and the other
  # column names are the vars from sample metadata that we use
  # na.or.complete removes rows with NAs, if no rows remain then correlation is NA
  corrResult <- data.table::as.data.table(cor(data1, data2, method = method, use='na.or.complete'), keep.rownames = T)

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

#' Correlation
#'
#' This function returns correlation coefficients for all columns in one data table against themselves.
#' 
#' @param data1 data.table with columns as variables. All columns must be numeric. One row per sample.
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return data.frame with correlation coefficients
setMethod("correlation", signature("data.table", "missing"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  # Check that all values are numeric
  if (!identical(veupathUtils::findNumericCols(data1), names(data1))) { stop("All columns in data1 must be numeric.")}

  ## Compute correlation
  # rownames and colnames should be the same in this case
  # na.or.complete removes rows with NAs, if no rows remain then correlation is NA
  corrResult <- cor(data1, method = method, use='na.or.complete')
  veupathUtils::logWithTime(paste0('Completed correlation with method=', method,'. Formatting results.'), verbose)

  ## Format results
  rowAndColNames <- expand.grid(rownames(corrResult), colnames(corrResult))
  deDupedRowAndColNames <- rowAndColNames[as.vector(upper.tri(corrResult)),]
  formattedCorrResult <- cbind(deDupedRowAndColNames, corrResult[upper.tri(corrResult)])
  colnames(formattedCorrResult) <- c("data1","data2","value")

  return(formattedCorrResult)
})

getDataMetadataType <- function(data) {
  if (inherits(data, 'AbundanceData')) {
    return('assay')
  } else if (inherits(data, 'SampleMetadata')) {
    return('sampleMetadata')
  } else {
    stop('data must be an AbundanceData or SampleMetadata object')
  }
}

## Helper function
# should this be s4?
buildCorrelationComputeResult <- function(corrResult, data1, data2 = NULL, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  method <- veupathUtils::matchArg(method)
  verbose <- veupathUtils::matchArg(verbose)
  if (!inherits(data1, c("AbundanceData", "SampleMetadata"))) {
    stop('data1 must be an AbundanceData or SampleMetadata object')
  }

  # both AbundanceData and SampleMetadata have these slots
  recordIdColumn <- data1@recordIdColumn
  ancestorIdColumns <- data1@ancestorIdColumns
  allIdColumns <- c(recordIdColumn, ancestorIdColumns)

  ## Format results
  # Construct the ComputeResult
  result <- new("ComputeResult")
  result@name <- 'correlation'
  result@recordIdColumn <- recordIdColumn
  result@ancestorIdColumns <- ancestorIdColumns
  statistics <- CorrelationResult(
    statistics = corrResult,
    data1Metadata = getDataMetadataType(data1),
    data2Metadata = ifelse(is.null(data2), getDataMetadataType(data1), getDataMetadataType(data2))
  )
  result@statistics <- statistics
  result@parameters <- paste0('method = ', method)

  # The resulting data should contain only the samples actually used.
  # this seems slightly complicated to generalize, and im not sure what we use it for anyhow
  #result@data <- abundances[, ..allIdColumns]
  #names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))

  validObject(result)
  veupathUtils::logWithTime(paste('Correlation computation completed with parameters recordIdColumn=', recordIdColumn, ', method = ', method), verbose)
  
  return(result)
}

#' Correlation of abundance data and metadata
#' 
#' This function returns the correlation of all columns of an AbundanceData object with appropriate columns of a SampleMetadata object.
#' 
#' @param data1 AbundanceData object. Will correlate abundance variables with specified variables in data2
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return a ComputeResult object
#' @export
setMethod("correlation", signature("AbundanceData", "missing"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  abundances <- getAbundances(data1, FALSE, FALSE)
  corrResult <- correlation(abundances, getSampleMetadata(data1, TRUE, FALSE), method, verbose)

  veupathUtils::logWithTime(paste("Received df table with", nrow(abundances), "samples and", (ncol(abundances)-1), "features with abundances."), verbose)

  result <- buildCorrelationComputeResult(corrResult, data1, data1@sampleMetadata, method, verbose)
  return(result)  
})

#' Self Correlation
#'
#' This function returns correlation coefficients for variables in one dataset against itself
#' 
#' @param data first dataset. An AbundanceData object or data.table
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setGeneric("selfCorrelation",
  function(data, method = c('spearman','pearson'), verbose = c(TRUE, FALSE), ...) standardGeneric("selfCorrelation"),
  signature = c("data")
)

#' Self Correlation
#'
#' This function returns correlation coefficients for variables in one AbundanceData object against itself
#' 
#' @param data An AbundanceData object
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setMethod("selfCorrelation", signature("AbundanceData"), function(data, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  abundances <- getAbundances(data, FALSE, FALSE)
  corrResult <- correlation(abundances, method=method, verbose=verbose)

  veupathUtils::logWithTime(paste("Received df table with", nrow(abundances), "samples and", (ncol(abundances)-1), "features with abundances."), verbose)

  result <- buildCorrelationComputeResult(corrResult, data, NULL, method, verbose)
  return(result)  
})

#' Self Correlation
#'
#' This function returns correlation coefficients for variables in one SampleMetadata object against itself
#' 
#' @param data SampleMetadata object
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setMethod("selfCorrelation", signature("SampleMetadata"), function(data, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  corrResult <- correlation(getSampleMetadata(data, TRUE, FALSE), method=method, verbose=verbose)

  veupathUtils::logWithTime(paste("Received df table with", nrow(data1), "samples and", (ncol(data)-1), "variables."), verbose)

  result <- buildCorrelationComputeResult(corrResult, data, NULL, method, verbose)
  return(result)
})

#' Self Correlation
#'
#' This function returns correlation coefficients for variables in one data.table against itself.
#' This is essentially an alias to the microbiomeComputations::correlation function.
#' 
#' @param data a data.table
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setMethod("selfCorrelation", signature("data.table"), function(data, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  correlation(data, method=method, verbose=verbose)
})

#' Correlation of two AbundanceData objects
#' 
#' This function returns the correlation of all columns of an AbundanceData object with appropriate columns of a second AbundanceData object.
#' 
#' @param data1 AbundanceData object. 
#' @param data2 AbundanceData object.
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @export
setMethod("correlation", signature("AbundanceData", "AbundanceData"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  abundances1 <- getAbundances(data1, FALSE, FALSE)
  abundances2 <- getAbundances(data2, FALSE, FALSE)
  corrResult <- correlation(abundances1, abundances2, method, verbose)

  veupathUtils::logWithTime(paste("Received first df table with", nrow(abundances1), "samples and", (ncol(abundances1)-1), "features with abundances."), verbose)
  veupathUtils::logWithTime(paste("Received second df table with", nrow(abundances2), "samples and", (ncol(abundances2)-1), "features with abundances."), verbose)

  result <- buildCorrelationComputeResult(corrResult, data1, data2, method, verbose)
  return(result) 
})






### DEPRECATED ###

#' Correlation of abundance data and metadata
#' 
#' This function returns the correlation of all columns of an AbundanceData object with appropriate columns of a SampleMetadata object.
#' This is deprecated. Please use correlation(AbundanceData) instead.
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
  veupathUtils::logWithTime(paste("Received df table with", nrow(df1), "samples and", (ncol(df1)-1), "features with abundances."), verbose)


  ## Compute correlation
  # Result has columns "data1", "data2", and correlationCoef
  corrResult <- correlation(df1[, -..allIdColumns], df2[, -..allIdColumns], method = method, verbose = verbose)


  ## Format results
  # Construct the ComputeResult
  result <- new("ComputeResult")
  result@name <- 'correlation'
  result@recordIdColumn <- recordIdColumn
  result@ancestorIdColumns <- ancestorIdColumns
  statistics <- CorrelationResult(
    statistics = corrResult,
    data1Metadata = "assay",
    data2Metadata = "sampleMetadata"
  )
  result@statistics <- statistics
  result@parameters <- paste0('method = ', method)


  # The resulting data should contain only the samples actually used.
  result@data <- df1[, ..allIdColumns]
  names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))


  validObject(result)
  veupathUtils::logWithTime(paste('Correlation computation completed with parameters recordIdColumn=', recordIdColumn, ', method = ', method), verbose)
  
  return(result)
})