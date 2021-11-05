#' Beta diversity
#'
#' This function returns pcoa coordinates calculated from the beta diversity dissimilarity matrix.
#' 
#' @param df data.frame with samples as rows, taxa as columns
#' @param sampleIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param method string defining the the beta diversity dissimilarity method. Accepted values are 'bray','jaccard', and 'jsd'
#' @param k integer determining the number of pcoa dimensions to return
#' @param verbose boolean indicating if timed logging is desired
#' @return who knows
#' @importFrom Rcpp sourceCpp
#' @importFrom vegan vegdist
#' @importFrom ape pcoa
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
betaDiv <- function(df,
                    sampleIdColumn,
                    method = c('bray','jaccard','jsd'),
                    k = 2,
                    verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", NROW(df), "samples and", (NCOL(df)-1), "taxa."), verbose)

    # Compute beta diversity using given dissimilarity method
    if (identical(method, 'bray') | identical(method, 'jaccard')) {

      dist <- vegan::vegdist(df[, -..sampleIdColumn], method=method, binary=TRUE)

    } else if (identical(method, 'jsd')) {

      dfMat <- matrix(as.numeric(unlist(df[, -..sampleIdColumn])), nrow=NROW(df))
      dist <- jsd(t(dfMat))
      dist <- as.dist(dist)

    } else {
      stop('Unaccepted dissimilarity method. Accepted methods are bray, jaccard, and jsd.')
    }
    veupathUtils::logWithTime("Computed dissimilarity matrix.", verbose)

    # Ordination
    ## Need to handle how this might err
    pcoa <- ape::pcoa(dist)
    dt <- data.table::as.data.table(pcoa$vectors)
    computeMessage <- paste("PCoA returned results for", NCOL(dt), "dimensions.")

    dt$SampleID <- df[[sampleIdColumn]]
    data.table::setcolorder(dt, c('SampleID'))
    data.table::setnames(dt,'SampleID', sampleIdColumn)
    veupathUtils::logWithTime("Finished ordination step.", verbose)

    # Extract percent variance
    eigenvecs <- pcoa$values$Relative_eig
    percentVar <- round(100*(eigenvecs / sum(eigenvecs)), 1)

    # Keep dims 1:k
    # We should keep the same number of percentVar values as cols in the data table. However, i think we're letting the user download lots of columns? So perhaps we shouldn't have k at all? A plot can use however many it needs.
    # For now returning data and percentVar for how much is in the plot.
    percentVar <- percentVar[1:k]
    keepCols <- c(sampleIdColumn,names(dt)[2:(k+1)])
    dt <- dt[, ..keepCols]

    #### Need to add back computed Variable Labels
    # Collect attributes
    entity <- veupathUtils::strSplit(sampleIdColumn,".", 4, 1)
    attr <- list('computationDetails' = computeMessage,
                 'parameters' = method,
                 'pcoaVariance' = percentVar,
                 'recordVariable' = sampleIdColumn)
    
    #### Make into a function? Need to get entity from variables and add display labels
    attr$computedVariableDetails <- list('id' = names(dt[, -..sampleIdColumn]),
                                         'entity' = entity,
                                         'displayLabel' = paste0(names(dt[, -..sampleIdColumn]), " ", percentVar, "%"),
                                         'isCollection' = FALSE)
    # Add entity to column names
    data.table::setnames(dt, names(dt[, -..sampleIdColumn]), paste0(entity,".",names(dt[, -..sampleIdColumn])))
    
    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    veupathUtils::logWithTime(paste('Beta diversity computation completed with parameters sampleIdColumn=', sampleIdColumn, ', method =', method, ', k =', k, ', verbose =', verbose), verbose)
    
    return(dt)
}

#' Beta diversity app
#'
#' This function returns the name of a json file with beta diversity results.
#' 
#' @param df data.frame with samples as rows, taxa as columns
#' @param sampleIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param methods vector of strings defining the the beta diversity dissimilarity methods to use. Must be a subset of c('bray','jaccard','jsd').
#' @param k integer determining the number of pcoa dimensions to return
#' @param verbose boolean indicating if timed logging is desired
#' @return we'll see.
#' @export
betaDivApp <- function(df,
                      sampleIdColumn,
                      methods = c('bray','jaccard','jsd'),
                      k = 2,
                      verbose = c(TRUE, FALSE)) {

    df <- data.table::setDT(df)

    verbose <- veupathUtils::matchArg(verbose)
    
    # Allow any number of methods to be inputted
    allMethods <- c('bray','jaccard','jsd')
    if (is.null(methods)) methods <- allMethods
    if (!all(methods %in% allMethods)) {
      stop("Unaccepted method found in 'methods' argument. 'methods' must be a subset of c('bray','jaccard','jsd').")
    }

    # Check that incoming df meets requirements - consider moving to a validateOTU function or similar
    if (!'data.table' %in% class(df)) {
      data.table::setDT(df)
    }
    if (!sampleIdColumn %in% names(df)) {
      stop("sampleIdColumn must exist as a column in df")
    }
    if (!all(unlist(lapply(df[, -..sampleIdColumn], is.numeric)))) {
      stop("All columns except the sampleIdColumn must be numeric")
    }

    appResults <- lapply(methods, betaDiv, df=df, sampleIdColumn=sampleIdColumn, k=k, verbose=verbose)
    

    # Write to json file - debating whether to keep this in here or move elsewhere. Makes testing easier
    # outFileName <- writeAppResultsToJson(appResults, 'betaDiv', verbose = verbose)

    return(appResults)

}