#' Ranked abundance
#'
#' This function returns abundances, ranked by a selected ranking function
#' 
#' @param df data.table with rows corresponding to records, one column indicating the record id, and all other columns corresponding to taxa abundance.
#' @param recordIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param method string defining the ranking strategy by which to order the taxa. Accepted values are 'median','max','q3',and 'variance'. Note that taxa that return a value of 0 for a given method will not be included in the results.
#' @param cutoff integer indicating the maximium number of taxa to be kept after ranking.
#' @param naToZero boolean indicating if NAs in numeric columns should be replaced with 0
#' @param verbose boolean indicating if timed logging is desired
#' @return data.table with columns recordIdColumn, followed by the top taxa as ranked by the given method, with no more taxa columns than specified by the cutoff argument.
#' @import veupathUtils
#' @import data.table
#' @export
rankedAbundance <- function(df, recordIdColumn, method = c('median','max','q3','variance'), cutoff=10, naToZero=c(TRUE, FALSE), verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    naToZero <- veupathUtils::matchArg(naToZero)
    verbose <- veupathUtils::matchArg(verbose)

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", nrow(df), "samples and", (ncol(df)-1), "taxa."), verbose)

    if (naToZero) {
      # Replace NA values with 0
      veupathUtils::setNaToZero(df)
      veupathUtils::logWithTime("Replaced NAs with 0", verbose)
    }

    # Reshape back to sample, taxonomicLevel, abundance
    formattedDT <- data.table::melt(df, measure.vars=colnames(df[, -..recordIdColumn]), variable.factor=F, variable.name='TaxonomicLevel', value.name="Abundance")

    rankedTaxa <- rankTaxa(formattedDT, method)

    # Extract top N taxa
    topN <- rankedTaxa[Abundance > 0, TaxonomicLevel]
    isCutoff <- FALSE
    if (length(topN) > cutoff) {
      topN <- topN[1:cutoff]
      computeMessage <- paste("Only returning top", cutoff, "results.")
      isCutoff <- TRUE
    }

    keepCols <- c(recordIdColumn, topN)
    dt = df[, ..keepCols]

    veupathUtils::logWithTime("Finished ranking taxa", verbose)
    
    # Collect attributes
    entity <- veupathUtils::strSplit(recordIdColumn,".", 4, 1)
    attr <- list('computationDetails' = computeMessage,
                 'parameters' = method,
                 'isCutoff' = isCutoff)
    
    #### Make into a function? Need to get entity from variables
    computedVariableDetails <- list('variableId' = unlist(lapply(names(dt[, -..recordIdColumn]), veupathUtils::strSplit, ".", 4, 2)),
                                    'entityId' = rep(entity, length(names(dt[, -..recordIdColumn]))),
                                    'dataType' = rep('NUMBER', length(names(dt[, -..recordIdColumn]))),
                                    'dataShape' = rep('CONTINUOUS', length(names(dt[, -..recordIdColumn]))),
                                    'isCollection' = TRUE)
    
    computedVariableMetadata <- list('displayRangeMin' = '0',
                                     'displayRangeMax' = '1',
                                     'collectionVariable' = list('collectionType' = 'abundance'))
    
    attr$computedVariable <- list('computedVariableDetails' = computedVariableDetails,
                                  'computedVariableMetadata' = computedVariableMetadata)
    
    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    veupathUtils::logWithTime(paste('Ranked abundance computation completed with parameters recordIdColumn=', recordIdColumn, ', method =', method, ', cutoff =', cutoff, ', naToZero = ', naToZero, ', verbose =', verbose), verbose)

    return(dt)
}

#' Ranked abundance app
#'
#' This function returns the name of a json file with ranked abundance results.
#' 
#' @param df data.frame with samples as rows, taxa as columns.
#' @param recordIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param methods vector of strings indicating ranking methods to use. Must be a subset of c('median','max','q3','variance').
#' @param cutoff integer indicating the maximium number of taxa to be kept after ranking.
#' @param naToZero boolean indicating if NAs in numeric columns should be replaced with 0
#' @param verbose boolean indicating if timed logging is desired.
#' @return name of a json file containing a list of data.tables, one for each method specified in methods. Each data.table contains a column for the recordIdColumn, top taxa columns, and an attribute "parameters" that records the method used.
#' @import veupathUtils
#' @export
rankedAbundanceApp <- function(df, recordIdColumn, methods=c('median','max','q3','variance'), cutoff=10, naToZero=c(TRUE, FALSE), verbose=c(TRUE, FALSE)) {

    naToZero <- veupathUtils::matchArg(naToZero)
    verbose <- veupathUtils::matchArg(verbose)

    allMethods <- c('median','max','q3','variance')
    # Allow any number of methods to be inputted
    if (is.null(methods)) methods <- allMethods
    if (!all(methods %in% allMethods)) {
      stop("Unaccepted method found in 'methods' argument. 'methods' must be a subset of c('median', 'max', 'q3','variance').")
    }
    
    # Check that incoming df meets requirements
    if (!'data.table' %in% class(df)) {
      data.table::setDT(df)
    }
    if (!recordIdColumn %in% names(df)) {
      stop("recordIdColumn must exist as a column in df")
    }
    if (!all(unlist(lapply(df[, -..recordIdColumn], is.numeric)))) {
      stop("All columns except the recordIdColumn must be numeric")
    }
    if (uniqueN(veupathUtils::strSplit(names(df), ".", ncol=2, index=1)) > 1) {
      stop("All entities must be identical")
    }

    appResults <- lapply(methods, rankedAbundance, df=df, recordIdColumn=recordIdColumn, cutoff=cutoff, naToZero=naToZero, verbose=verbose)

    # Write to json file
    # outFileName <- writeAppResultsToJson(appResults, recordIdColumn = recordIdColumn, 'rankedAbundance', verbose = verbose)

    return(appResults)

}