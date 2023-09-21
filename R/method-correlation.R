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
#' @param data2 SampleMetadata object
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @param variables VariableMetadataList object with variables corresponding to columns in the SampleMetadata object.
#' @param data1Name string name for column in the correlation results that corresponds to variables from data1
#' @param data2Name string name for column in the correlation results that corresponds to variables from data2
#' @export
setMethod("correlation", signature("AbundanceData", "SampleMetadata"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE), variables = NULL, data1Name = "data1", data2Name = "data2") {
    df1 <- data1@data
    recordIdColumn1 <- data1@recordIdColumn
    naToZero <- data1@imputeZero
    ancestorIdColumns1 <- data1@ancestorIdColumns
    allIdColumns1 <- c(recordIdColumn1, ancestorIdColumns1)
    df2 <- data2@data
    recordIdColumn2 <- data2@recordIdColumn
    ancestorIdColumns2 <- data2@ancestorIdColumns
    allIdColumns2 <- c(recordIdColumn2, ancestorIdColumns2)


    ## Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    # Check that incoming abundance data meets requirements
    if (!'data.table' %in% class(df1)) {
      data.table::setDT(df1)
    }

    # Check that incoming metadata meets requirements
    if (!'data.table' %in% class(df2)) {
      data.table::setDT(df2)
    }

    # Check that variables is the correct class
    if (!is.null(variables) && !'VariableMetadataList' %in% class(variables)) {
      stop("The variables argument must be a VariableMetadataList if defined.")
    }

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", nrow(df1), "samples and", (ncol(df1)-1), "taxa."), verbose)

    # Replace NA values with 0
    if (naToZero) {
      veupathUtils::setNaToZero(df1)
      veupathUtils::logWithTime("Replaced NAs with 0.", verbose)
    }

    # Ensure the sample IDs, and their order, are the same in both datasets
    if (!identical(df1[[recordIdColumn1]], df2[[recordIdColumn2]])) {
      if (!identical(union(df1[[recordIdColumn1]], df2[[recordIdColumn2]]), intersect(df1[[recordIdColumn1]], df2[[recordIdColumn2]]))) {
        stop("Sample IDs do not match between data1 and data2.")
      }
      # Reorder df2 (sample metadata) so that it matches the order in df
      sampleMetadata <- df2[match(df1[[recordIdColumn1]], df2[[recordIdColumn2]]), ]
      veupathUtils::logWithTime("Reordered sampleMetadata rows based on abundance data sample id order.", verbose)
    }

    ## Extract appropriate metadata columns
    if (!is.null(variables)) {
      # Use the data shape to find appropriate columns in the metadata
      inputMetadataCols <- veupathUtils::findColNamesByPredicate(variables, function(x) {identical(x@dataShape@value, "CONTINUOUS")})

      # If any are dates, we'll need to coerce them to numeric for the correlation
      dateMetadataCols <- veupathUtils::findColNamesByPredicate(variables, function(x) {identical(x@dataShape@value, "CONTINUOUS") && identical(x@dataType@value, "DATE")})
      if (!!length(dateMetadataCols)) {
        veupathUtils::logWithTime("Converting date variables to numeric", verbose)
        df2[, (dateMetadataCols) := lapply(.SD, as.numeric), .SDcols = dateMetadataCols]
      }

    } else {
      warning("No variable metadata supplied. Using all numeric columns of sampleMetadata for correlation.")
      inputMetadataCols <- veupathUtils::findNumericCols(df2[, -..allIdColumns2])
    }
    if (length(inputMetadataCols) < 1) {
      stop("correlation requires at least one continuous metadata variable.")
    }


    ## Compute correlation
    # Result has columns "data1", "data2", and correlationCoef
    corrResult <- correlation(df1[, -..allIdColumns1], df2[, ..inputMetadataCols], method = method, verbose = verbose)


    ## Format results
    # Rename the correlation results columns so that they are a little more friendly
    names(corrResult)[names(corrResult) == "data1"] <- data1Name
    names(corrResult)[names(corrResult) == "data2"] <- data2Name

    
    # Construct the ComputeResult
    result <- new("ComputeResult")
    result@name <- 'correlation'
    result@recordIdColumn <- recordIdColumn1
    result@ancestorIdColumns <- ancestorIdColumns1
    result@statistics <- corrResult
    result@parameters <- paste0('method = ', method)


    # The resulting data should contain only the samples actually used.
    result@data <- df1[, ..allIdColumns1]
    names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))


    validObject(result)
    veupathUtils::logWithTime(paste('Correlation computation completed with parameters recordIdColumn=', recordIdColumn1, ', method = ', method), verbose)
    
    return(result)
})