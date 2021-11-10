#' Beta diversity
#'
#' This function returns pcoa coordinates calculated from the beta diversity dissimilarity matrix.
#' 
#' @param df data.table with rows corresponding to records, one column indicating the record id, and all other columns corresponding to taxa abundance.
#' @param recordIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param method string defining the the beta diversity dissimilarity method. Accepted values are 'bray','jaccard', and 'jsd'
#' @param k integer determining the number of pcoa axes to return
#' @param verbose boolean indicating if timed logging is desired
#' @return data.table with a column for the recordIdColumn, as well as columns corresponding to pcoa axes.
#' @importFrom Rcpp sourceCpp
#' @importFrom vegan vegdist
#' @importFrom ape pcoa
#' @importFrom stringi stri_replace_all_fixed
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
betaDiv <- function(df,
                    recordIdColumn,
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

      dist <- try(vegan::vegdist(df[, -..recordIdColumn], method=method, binary=TRUE))

    } else if (identical(method, 'jsd')) {

      dfMat <- matrix(as.numeric(unlist(df[, -..recordIdColumn])), nrow=NROW(df))
      dist <- try({dist <- jsd(t(dfMat)); dist <- as.dist(dist)})

    } else {
      stop('Unaccepted dissimilarity method. Accepted methods are bray, jaccard, and jsd.')
    }
    
    # Handle errors or return positive computeMessage
    if (veupathUtils::is.error(dist)) {
      
      computeMessage <- paste('Error: beta diversity', method, 'failed:', attr(dist,'condition')$message)
      
      # Return only recordIdColumn and expected attributes
      dt <- df[, ..recordIdColumn]
      
      attr <- list('computationDetails' = computeMessage,
                   'parameters' = character(),
                   'recordVariable' = character(),
                   'pcoaVariance' = numeric())
      attr$computedVariableDetails <- list('variableId' = character(),
                                           'entityId' = character(),
                                           'dataType' = character(),
                                           'dataShape' = character(),
                                           'displayLabel' = character(),
                                           'isCollection' = logical())
      veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)
      veupathUtils::logWithTime(paste('Beta diversity computation FAILED with parameters recordIdColumn=', recordIdColumn, ', method=', method, ', k=', k , ', verbose =', verbose), verbose)
      
      return(dt)
      
    } else {
      veupathUtils::logWithTime("Computed dissimilarity matrix.", verbose)
    }

    # Ordination
    pcoa <- ape::pcoa(dist)
    dt <- data.table::as.data.table(pcoa$vectors)
    # Remove dots from names
    data.table::setnames(dt, stringi::stri_replace_all_fixed(names(dt),".",""))
    computeMessage <- paste("PCoA returned results for", NCOL(dt), "dimensions.")

    dt[[recordIdColumn]] <- df[[recordIdColumn]]
    data.table::setcolorder(dt, recordIdColumn)
    veupathUtils::logWithTime("Finished ordination step.", verbose)

    # Extract percent variance
    eigenvecs <- pcoa$values$Relative_eig
    percentVar <- round(100*(eigenvecs / sum(eigenvecs)), 1)

    # Keep dims 1:k
    # We should keep the same number of percentVar values as cols in the data table. However, i think we're letting the user download lots of columns? So perhaps we shouldn't have k at all? A plot can use however many it needs.
    # For now returning data and percentVar for how much is in the plot.
    percentVar <- percentVar[1:k]
    keepCols <- c(recordIdColumn,names(dt)[2:(k+1)])
    dt <- dt[, ..keepCols]

    #### Need to add back computed Variable Labels
    # Collect attributes
    entity <- veupathUtils::strSplit(recordIdColumn,".", 4, 1)
    attr <- list('computationDetails' = computeMessage,
                 'parameters' = method,
                 'pcoaVariance' = percentVar,
                 'recordVariable' = recordIdColumn)
    
    #### Make into a function? Need to get entity from variables and add display labels
    attr$computedVariableDetails <- list('variableId' = names(dt[, -..recordIdColumn]),
                                         'entityId' = rep(entity, length(names(dt[, -..recordIdColumn]))),
                                         'dataType' = rep('NUMBER', length(names(dt[, -..recordIdColumn]))),
                                         'dataShape' = rep('CONTINUOUS', length(names(dt[, -..recordIdColumn]))),
                                         'displayLabel' = paste0(names(dt[, -..recordIdColumn]), " ", sprintf(percentVar,fmt = '%#.1f'), "%"),
                                         'isCollection' = FALSE)
    # Add entity to column names
    data.table::setnames(dt, names(dt[, -..recordIdColumn]), paste0(entity,".",names(dt[, -..recordIdColumn])))
    
    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    veupathUtils::logWithTime(paste('Beta diversity computation completed with parameters recordIdColumn=', recordIdColumn, ', method =', method, ', k =', k, ', verbose =', verbose), verbose)
    
    return(dt)
}

#' Beta diversity app
#'
#' This function returns the name of a json file with beta diversity results.
#' 
#' @param df data.frame with samples as rows, taxa as columns
#' @param recordIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param methods vector of strings defining the the beta diversity dissimilarity methods to use. Must be a subset of c('bray','jaccard','jsd').
#' @param k integer determining the number of pcoa dimensions to return
#' @param verbose boolean indicating if timed logging is desired
#' @return name of a json file containing a list of data.tables, one for each method specified in methods. Each data.table contains a column for the recordIdColumn, as well as columns corresponding to pcoa axes.
#' @export
betaDivApp <- function(df,
                      recordIdColumn,
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
    if (!recordIdColumn %in% names(df)) {
      stop("recordIdColumn must exist as a column in df")
    }
    if (!all(unlist(lapply(df[, -..recordIdColumn], is.numeric)))) {
      stop("All columns except the recordIdColumn must be numeric")
    }
    if (uniqueN(veupathUtils::strSplit(names(df), ".", ncol=2, index=1)) > 1) {
      stop("All entities must be identical")
    }

    appResults <- lapply(methods, betaDiv, df=df, recordIdColumn=recordIdColumn, k=k, verbose=verbose)
    

    # Write to json file - debating whether to keep this in here or move elsewhere. Makes testing easier
    # outFileName <- writeAppResultsToJson(appResults, 'betaDiv', verbose = verbose)

    return(appResults)

}