#' Predicate Factory
#' 
#' This function creates a predicate function based on a string defining the type of predicate to run and a numeric value. 
#' The currently supported types are 'proportionNonZero', 'variance' and 'sd'. The numeric value associated
#' with each predicate type is the threshold for the predicate to be true.
#' 
#' @param predicateType string defining the type of predicate to run. The currently supported values are 'proportionNonZero', 'variance' and 'sd'
#' @param threshold numeric value associated with the predicate type
#' @return Function returning a boolean indicating if a feature should be included (TRUE) or excluded (FALSE)
#' @export
setGeneric("predicateFactory",
  function(predicateType, threshold) standardGeneric("predicateFactory"),
  signature = c("predicateType", "threshold")
)


#' Predicate Factory
#'
#' This function creates a predicate function based on a string defining the type of predicate to run and a numeric value. 
#' The currently supported types are 'proportionNonZero', 'variance' and 'sd'. The numeric value associated
#' with each predicate type is the threshold for the predicate to be true.
#' 
#' @param predicateType string defining the type of predicate to run. The currently supported values are 'proportionNonZero', 'variance' and 'sd'
#' @param threshold numeric value associated with the predicate type
#' @return Function returning a boolean indicating if a feature should be included (TRUE) or excluded (FALSE)
#' @export
setMethod("predicateFactory", signature("character", "numeric"), function(predicateType = c('proportionNonZero', 'variance', 'sd'), threshold = 0.5) {
  predicateType <- veupathUtils::matchArg(predicateType)

  if (predicateType == 'proportionNonZero') {
    if (threshold < 0 | threshold > 1) {
      stop('threshold must be between 0 and 1 for proportionNonZero')
    }
    return(function(x){sum(x > 0) >= length(x) * threshold})

  } else if (predicateType == 'variance') {
    if (threshold < 0) {
      stop('threshold must be greater than 0 for variance')
    }
    return(function(x){var(x) > threshold})

  } else if (predicateType == 'sd') {
    if (threshold < 0) {
      stop('threshold must be greater than 0 for sd')
    }
    return(function(x){sd(x) > threshold})
  }

})

#' Correlation
#'
#' This function returns correlation coefficients for variables in one dataset against variables in a second dataset
#' 
#' @param data1 first dataset. An AbundanceData object or data.table
#' @param data2 second dataset. A SampleMetadata object (if data1 is class AbundanceData) or a data.table
#' @param method string defining the type of correlation to run. The currently supported values are specific to the class of data1 and data2.
#' @param verbose boolean indicating if timed logging is desired
#' @return data.frame with correlation coefficients or a ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
setGeneric("correlation",
  function(data1, data2, method, verbose = c(TRUE, FALSE), ...) standardGeneric("correlation"),
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
#' @importFrom Hmisc rcorr
#' @return data.frame with correlation coefficients
setMethod("correlation", signature("data.table", "data.table"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  method <- veupathUtils::matchArg(method)
  verbose <- veupathUtils::matchArg(verbose)

  # Check that the number of rows match.
  if (!identical(nrow(data1), nrow(data2))) {
    stop("data1 and data2 must have the same number of rows.")
  }

  # Check that all values are numeric
  if (!identical(veupathUtils::findNumericCols(data1), names(data1))) { stop("All columns in data1 must be numeric.")}
  if (!identical(veupathUtils::findNumericCols(data2), names(data2))) { stop("All columns in data2 must be numeric.")}


  ## Compute correlation
  #corrResult <- data.table::as.data.table(cor(data1, data2, method = method, use='na.or.complete'), keep.rownames = T)
  lastData1ColIndex <- length(data1)
  firstData2ColIndex <- length(data1) + 1
  corrResult <- Hmisc::rcorr(as.matrix(data1), as.matrix(data2), type = method)
  # this bc Hmisc::rcorr cbinds the two data.tables and runs the correlation
  # so we need to extract only the relevant values
  # of note, seems subsetting a matrix to have a single row/ column drops rownames by default. adding drop=F to prevent this.
  pVals <- data.table::as.data.table(corrResult$P[1:lastData1ColIndex, firstData2ColIndex:length(colnames(corrResult$P)), drop = F], keep.rownames = T)
  corrResult <- data.table::as.data.table(corrResult$r[1:lastData1ColIndex, firstData2ColIndex:length(colnames(corrResult$r)), drop = F], keep.rownames = T)

  veupathUtils::logWithTime(paste0('Completed correlation with method=', method,'. Formatting results.'), verbose)


  ## Format results
  meltedCorrResult <- melt(corrResult, id.vars=c('rn'))
  meltedPVals <- melt(pVals, id.vars=c('rn'))
  formattedCorrResult <- data.frame(
    data1 = meltedCorrResult[['rn']],
    data2 = meltedCorrResult[['variable']],
    correlationCoef = meltedCorrResult[['value']],
    # should we do a merge just to be sure?
    pValue = meltedPVals[['value']]
  )

  return(formattedCorrResult)
})

#' Correlation
#'
#' This function returns correlation coefficients for all columns in one data table against themselves.
#' 
#' @param data1 data.table with columns as variables. All columns must be numeric. One row per sample.
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman', 'pearson' and 'sparcc'
#' @param verbose boolean indicating if timed logging is desired
#' @return data.frame with correlation coefficients
#' @importFrom Hmisc rcorr
#' @importFrom SpiecEasi pval.sparccboot
#' @importFrom SpiecEasi sparccboot
#' @importFrom SpiecEasi sparcc
setMethod("correlation", signature("data.table", "missing"), function(data1, data2, method = c('spearman','pearson','sparcc'), verbose = c(TRUE, FALSE)) {
  method <- veupathUtils::matchArg(method)
  verbose <- veupathUtils::matchArg(verbose)

  # Check that all values are numeric
  if (!identical(veupathUtils::findNumericCols(data1), names(data1))) { stop("All columns in data1 must be numeric.")}

  ## Compute correlation
  # rownames and colnames should be the same in this case
  # keep matrix for now so we can use lower.tri later, expand.grid will give us the needed data.frame
  if (method == 'sparcc') {

    # this is a local alias for the sparcc function from the SpiecEasi namespace that do.call can find
    # if we need to customize how we call sparcc in the future, we'll need to do something like the following
    # except the empty list would have args for the sparcc function
    #sparcc <- get("sparcc", asNamespace("SpiecEasi"))
    #statisticperm=function(data, indices) do.call("sparcc", c(list(apply(data[indices,], 2, sample)), list()))$Cor
    #statisticboot=function(data, indices) do.call("sparcc", c(list(data[indices,,drop=FALSE]), list()))$Cor

    # sub-sampled `statistic` functions to pass to `boot::boot` that we can use to find pvalues
    # statisticboot = function which takes data and bootstrap sample indices of the bootstapped correlation matrix
    # statisticperm = function which takes data and permutated sample indices of the null correlation matrix
    # the SpiecEasi defaults for these functions return the upper triangle of the correlation matrix. We do that manually later instead.
    statisticperm = function(data, indices) SpiecEasi::sparcc(apply(data[indices,], 2, sample))$Cor
    statisticboot = function(data, indices) SpiecEasi::sparcc(data[indices,,drop=FALSE])$Cor

    # calling the bootstrap version of sparcc and finding pvalues
    result <- SpiecEasi::pval.sparccboot(SpiecEasi::sparccboot(data1, statisticboot = statisticboot, statisticperm = statisticperm, R = 100))
    # making sure results are formatted correctly for downstream use
    pVals <- matrix(result$pvals, nrow = ncol(data1), ncol = ncol(data1), byrow = T)
    corrResult <- result$cors
    rownames(corrResult) <- colnames(corrResult) <- colnames(data1)
  
  } else {
  
    corrResult <- Hmisc::rcorr(as.matrix(data1), type = method)
    pVals <- corrResult$P
    corrResult <- corrResult$r
  
  }

  veupathUtils::logWithTime(paste0('Completed correlation with method=', method,'. Formatting results.'), verbose)

  ## Format results
  rowAndColNames <- expand.grid(rownames(corrResult), colnames(corrResult))
  deDupedRowAndColNames <- rowAndColNames[as.vector(upper.tri(corrResult)),]
  formattedCorrResult <- cbind(deDupedRowAndColNames, corrResult[upper.tri(corrResult)])
  formattedCorrResult <- cbind(formattedCorrResult, pVals[upper.tri(pVals)])
  colnames(formattedCorrResult) <- c("data1","data2","correlationCoef","pValue")

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
buildCorrelationComputeResult <- function(corrResult, data1, data2 = NULL, method = c('spearman','pearson','sparcc'), verbose = c(TRUE, FALSE)) {
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
#' @param proportionNonZeroThreshold numeric threshold to filter features by proportion of non-zero values across samples
#' @param varianceThreshold numeric threshold to filter features by variance across samples
#' @param stdDevThreshold numeric threshold to filter features by standard deviation across samples
#' @return a ComputeResult object
#' @export
setMethod("correlation", signature("AbundanceData", "missing"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE), proportionNonZeroThreshold = 0.5, varianceThreshold = 0, stdDevThreshold = 0) {
  #prefilters applied
  data1 <- pruneFeatures(data1, predicateFactory('proportionNonZero', proportionNonZeroThreshold), verbose)
  data1 <- pruneFeatures(data1, predicateFactory('variance', varianceThreshold), verbose)
  data1 <- pruneFeatures(data1, predicateFactory('sd', stdDevThreshold), verbose)
  
  abundances <- getAbundances(data1, FALSE, FALSE, verbose)
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
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman','pearson' and 'sparcc'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setGeneric("selfCorrelation",
  function(data, method = c('spearman','pearson','sparcc'), verbose = c(TRUE, FALSE), ...) standardGeneric("selfCorrelation"),
  signature = c("data")
)

#' Self Correlation
#'
#' This function returns correlation coefficients for variables in one AbundanceData object against itself
#' 
#' @param data An AbundanceData object
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman','pearson' and 'sparcc'
#' @param verbose boolean indicating if timed logging is desired
#' @param proportionNonZeroThreshold numeric threshold to filter features by proportion of non-zero values across samples
#' @param varianceThreshold numeric threshold to filter features by variance across samples
#' @param stdDevThreshold numeric threshold to filter features by standard deviation across samples
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setMethod("selfCorrelation", signature("AbundanceData"), function(data, method = c('spearman','pearson','sparcc'), verbose = c(TRUE, FALSE), proportionNonZeroThreshold = 0.5, varianceThreshold = 0, stdDevThreshold = 0) {
  #prefilters applied
  data <- pruneFeatures(data, predicateFactory('proportionNonZero', proportionNonZeroThreshold), verbose)
  data <- pruneFeatures(data, predicateFactory('variance', varianceThreshold), verbose)
  data <- pruneFeatures(data, predicateFactory('sd', stdDevThreshold), verbose)

  abundances <- getAbundances(data, FALSE, FALSE, verbose)
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
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman', 'pearson' and 'sparcc'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setMethod("selfCorrelation", signature("SampleMetadata"), function(data, method = c('spearman','pearson','sparcc'), verbose = c(TRUE, FALSE)) {
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
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman', 'pearson' and 'sparcc'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @export
setMethod("selfCorrelation", signature("data.table"), function(data, method = c('spearman','pearson','sparcc'), verbose = c(TRUE, FALSE)) {
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
#' @param proportionNonZeroThreshold numeric threshold to filter features by proportion of non-zero values across samples
#' @param varianceThreshold numeric threshold to filter features by variance across samples
#' @param stdDevThreshold numeric threshold to filter features by standard deviation across samples
#' @return ComputeResult object
#' @export
setMethod("correlation", signature("AbundanceData", "AbundanceData"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE), proportionNonZeroThreshold = 0.5, varianceThreshold = 0, stdDevThreshold = 0) {
  #prefilters applied
  data1 <- pruneFeatures(data1, predicateFactory('proportionNonZero', proportionNonZeroThreshold), verbose)
  data1 <- pruneFeatures(data1, predicateFactory('variance', varianceThreshold), verbose)
  data1 <- pruneFeatures(data1, predicateFactory('sd', stdDevThreshold), verbose)
  data2 <- pruneFeatures(data2, predicateFactory('proportionNonZero', proportionNonZeroThreshold), verbose)
  data2 <- pruneFeatures(data2, predicateFactory('variance', varianceThreshold), verbose)
  data2 <- pruneFeatures(data2, predicateFactory('sd', stdDevThreshold), verbose)
  
  abundances1 <- getAbundances(data1, FALSE, TRUE, verbose)
  abundances2 <- getAbundances(data2, FALSE, TRUE, verbose)

  # empty samples removed from data by getAbundances, means we need to keep samples common to both datasets and remove ids
  # get id col names
  recordIdColumn <- data1@recordIdColumn
  allIdColumns <- c(recordIdColumn, data1@ancestorIdColumns)
  # should we verify that ids are the same in both datasets?

  # remove samples that are not common
  commonSamples <- intersect(abundances1[[recordIdColumn]], abundances2[[recordIdColumn]])
  abundances1 <- abundances1[abundances1[[recordIdColumn]] %in% commonSamples, ]
  abundances2 <- abundances2[abundances2[[recordIdColumn]] %in% commonSamples, ]

  # remove ids
  abundances1 <- abundances1[, -..allIdColumns]
  abundances2 <- abundances2[, -..allIdColumns]  

  corrResult <- correlation(abundances1, abundances2, method, verbose)

  veupathUtils::logWithTime(paste("Received first df table with", nrow(abundances1), "samples and", (ncol(abundances1)-1), "features with abundances."), verbose)
  veupathUtils::logWithTime(paste("Received second df table with", nrow(abundances2), "samples and", (ncol(abundances2)-1), "features with abundances."), verbose)

  result <- buildCorrelationComputeResult(corrResult, data1, data2, method, verbose)
  return(result) 
})
