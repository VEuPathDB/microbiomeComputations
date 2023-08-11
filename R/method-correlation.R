#' Correlation
#'
#' This function returns correlation coefficients and significance values for all-to-all correlation of variables in two datasets.
#' 
#' @param data whatAmI? Abundance data i guess
#' @param secondMatrix data table or something (assay 2) to correlate against abundace data. Optional
#' @param method string. Can be X, Y, or Z
#' @param verbose boolean
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
setGeneric("correlation",
  function(data, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) standardGeneric("correlation"),
  signature = c("data")
)

#'@export
setMethod("correlation", signature("AbundanceData"), function(data, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
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

    ## TODO Check that samples match. If not, reorder so they do

    ## Find numeric metadata
    numericMetadataCols <- veupathUtils::findNumericCols(sampleMetadata)

    ## (remember here samples are rows and vars are cols). Abundance data should be same!

    ## Compute correlation
    # Resulting data table has column "rn" = row names of the correlation matrix (so taxa names), and the other
    # column names are the vars from sample metadata that we use
    corrResult <- data.table::as.data.table(cor(df[, -..allIdColumns], sampleMetadata[, ..numericMetadataCols], method = method), keep.rownames = T)

    melted <- melt(corrResult, id.vars=c('rn'))
    # Output: rn = taxa or gene name, variable = metadata variable, value = correlation coefficient
    data.table::setnames(melted, c('rn','variable','value'), c('var1','var2','corrCoeff'))

    # Format results
    statistics <- data.frame(
      var1 = melted[['var1']],
      var2 = melted[['var2']],
      corrCoeff = melted[['corrCoeff']]
      # pValue = TODO,
      # adjustedPValue = TODO,
    )

    veupathUtils::logWithTime(paste0('Completed method=',method,'. Formatting results.'), verbose)
    
    ## Construct the ComputeResult
    result <- new("ComputeResult")
    result@name <- 'correlation'
    result@recordIdColumn <- recordIdColumn
    result@ancestorIdColumns <- ancestorIdColumns
    result@statistics <- statistics
    result@parameters <- paste0('method = ', method)
    # result@droppedColumns <- droppedColumns   TODO - should be added when we handle NAs


    # The resulting data should contain only the samples actually used.
    result@data <- df[, ..allIdColumns]
    names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))


    validObject(result)
    veupathUtils::logWithTime(paste('Correlation computation completed with parameters recordIdColumn=', recordIdColumn, ', method = ', method), verbose)
    
    return(result)
})