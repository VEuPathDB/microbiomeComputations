#' Correlation
#'
#' This function returns correlation coefficients for all abundance variables vs all continuous sample metadata variables.
#' Soon this function will be expanded to take multiple abundance datasets, and to run all-to-all correlations
#' for abundance variables.
#' 
#' @param  AbundanceData object
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman' and 'pearson'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
setGeneric("correlation",
  function(data1, data2, method = c('spearman','pearson'), ..., verbose = c(TRUE, FALSE)) standardGeneric("correlation"),
  signature = c("data1","data2")
)

setMethod("correlation", signature("data.table", "data.table"), function(data1, data2, method = c('spearman','pearson'), verbose = c(TRUE, FALSE)) {
  # data table correlation 

  # Check that the number of rows match. We'll do all cols by all cols

  # Check that all values are numeric

  ## Compute correlation
  # Resulting data table has column "rn" = row names of the correlation matrix (so taxa names), and the other
  # column names are the vars from sample metadata that we use
  corrResult <- data.table::as.data.table(cor(data1, data2, method = method), keep.rownames = T)

  veupathUtils::logWithTime(paste0('Completed correlation with method=',method,'. Formatting results.'), verbose)


  ## Format results
  meltedCorrResult <- melt(corrResult, id.vars=c('rn'))
  formattedCorrResult <- data.frame(
    data1 = meltedCorrResult[['rn']],
    data2 = meltedCorrResult[['variable']],
    correlationCoef = meltedCorrResult[['value']]
  )

  return(formattedCorrResult)

})

#' This is where the names should go, because they know what's data1 and data 2 and can do the appropriate renaming
#'@export
setMethod("correlation", signature("AbundanceData", "SampleMetadata"), function(data1, data2, method = c('spearman','pearson'), variables = NULL, data2Name = NULL,data1Name = NULL, verbose = c(TRUE, FALSE)) {
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

    # Check that incoming df meets requirements
    if (!'data.table' %in% class(df1)) {
      data.table::setDT(df1)
    }

    # Check that incoming metadata meets requirements
    if (!'data.table' %in% class(df2)) {
      data.table::setDT(df2)
    }

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", nrow(df1), "samples and", (ncol(df1)-1), "taxa."), verbose)

    # Replace NA values with 0
    if (naToZero) {
      veupathUtils::setNaToZero(df1)
      veupathUtils::logWithTime("Replaced NAs with 0", verbose)
    }

    # Subset to only use intersect of sample ids or just err

    # Check that the order of sample ids match. If not, reorder so they do
    if (!identical(df1[[recordIdColumn1]], df2[[recordIdColumn2]])) {
      # Reorder df2 (sample metadata) so that it matches the order in df
      sampleMetadata <- df2[match(df1[[recordIdColumn1]], df2[[recordIdColumn2]]), ]
      veupathUtils::logWithTime("Reordered sampleMetadata rows based on abundance data sample id order.", verbose)
    }

    ## Extract appropriate metadata columns
    if (!is.null(variables)) {
      # Use the data shape to find appropriate columns in the metadata
      inputMetadataCols <- veupathUtils::findColNamesByPredicate(variables, function(x) {identical(x@dataShape@value, "CONTINUOUS")})

      # If any are dates, we'll need to coerce them to numeric for the correlation, and then back again.
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
    # Resulting data table has column "rn" = row names of the correlation matrix (so taxa names), and the other
    # column names are the vars from sample metadata that we use
    corrResult <- correlation(df1[, -..allIdColumns1], df2[, ..inputMetadataCols], method = method, verbose = verbose)


    ## Format results
    print(class(data1)[1])
    print(names(corrResult))
    print(data1Name)
    # names(corrResult)[names(corrResult) == "data1"] <- ifelse(!!length(data1Name), data1Name, class(data1)[1])
    names(corrResult)[names(corrResult) == "data1"] <- ifelse(!is.null(data1Name), "test", "wrong")
    names(corrResult)[names(corrResult) == "data2"] <- ifelse(!is.null(data2Name), data2Name, class(data2)[1])
    print(head(corrResult))

    
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