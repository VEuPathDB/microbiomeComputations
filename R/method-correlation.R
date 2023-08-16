#' Correlation
#'
#' This function returns correlation coefficients for all abundance variables vs all continuous sample metadata variables.
#' Soon this function will be expanded to take multiple abundance datasets, and to run all-to-all correlations
#' for abundance variables.
#' 
#' @param data AbundanceData object
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
setGeneric("correlation",
  function(data, method = c('spearman','pearson'), variables = NULL, verbose = c(TRUE, FALSE)) standardGeneric("correlation"),
  signature = c("data")
)

#'@export
setMethod("correlation", signature("AbundanceData"), function(data, method = c('spearman','pearson'), variables = NULL, verbose = c(TRUE, FALSE)) {
    df <- data@data
    recordIdColumn <- data@recordIdColumn
    naToZero <- data@imputeZero
    ancestorIdColumns <- data@ancestorIdColumns
    allIdColumns <- c(recordIdColumn, ancestorIdColumns)
    sampleMetadata <- data@sampleMetadata

    ## Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    # Check that incoming df meets requirements
    if (!'data.table' %in% class(df)) {
      data.table::setDT(df)
    }

    # Check that incoming metadata meets requirements
    if (!'data.table' %in% class(sampleMetadata)) {
      data.table::setDT(sampleMetadata)
    }

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", nrow(df), "samples and", (ncol(df)-1), "taxa."), verbose)

    # Replace NA values with 0
    if (naToZero) {
      veupathUtils::setNaToZero(df)
      veupathUtils::logWithTime("Replaced NAs with 0", verbose)
    }

    # Check that the order of sample ids match. If not, reorder so they do
    if (!identical(df[[recordIdColumn]], sampleMetadata[[recordIdColumn]])) {
      # Reorder sampleMetadata so that it matches the order in df
      sampleMetadata <- sampleMetadata[match(df[[recordIdColumn]],sampleMetadata[[recordIdColumn]]), ]
      veupathUtils::logWithTime("Reordered sampleMetadata rows based on abundance data sample id order.", verbose)
    }

    ## Extract appropriate metadata columns
    if (!is.null(variables)) {
      # Use the data shape to find appropriate columns in the metadata
      inputMetadataCols <- veupathUtils::findColNamesByPredicate(variables, function(x) {identical(x@dataShape@value, "CONTINUOUS")})

    } else {
      warning("No variable metadata supplied. Using all numeric columns of sampleMetadata for correlation.")
      inputMetadataCols <- veupathUtils::findNumericCols(sampleMetadata)
    }
    if (length(inputMetadataCols) < 1) {
      stop("correlation requires at least one continuous metadata variable.")
    }


    ## Compute correlation
    # Resulting data table has column "rn" = row names of the correlation matrix (so taxa names), and the other
    # column names are the vars from sample metadata that we use
    corrResult <- data.table::as.data.table(cor(df[, -..allIdColumns], sampleMetadata[, ..inputMetadataCols], method = method), keep.rownames = T)

    veupathUtils::logWithTime(paste0('Completed method=',method,'. Formatting results.'), verbose)


    ## Format results
    meltedCorrResult <- melt(corrResult, id.vars=c('rn'))
    statistics <- data.frame(
      var1 = meltedCorrResult[['rn']],
      var2 = meltedCorrResult[['variable']],
      correlationCoef = meltedCorrResult[['value']]
    )
    
    # Construct the ComputeResult
    result <- new("ComputeResult")
    result@name <- 'correlation'
    result@recordIdColumn <- recordIdColumn
    result@ancestorIdColumns <- ancestorIdColumns
    result@statistics <- statistics
    result@parameters <- paste0('method = ', method)


    # The resulting data should contain only the samples actually used.
    result@data <- df[, ..allIdColumns]
    names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))


    validObject(result)
    veupathUtils::logWithTime(paste('Correlation computation completed with parameters recordIdColumn=', recordIdColumn, ', method = ', method), verbose)
    
    return(result)
})