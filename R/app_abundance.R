#' Ranked abundance
#'
#' This function returns abundances, ranked by a selected ranking function
#' 
#' @param otu data.frame with samples as rows, taxa as columns
#' @param method string defining the ranking strategy by which to order the taxa. Accepted values are 'median','max','q3',and 'var'.
#' @param cutoff integer indicating the maximium number of taxa to be kept after ranking.
#' @param verbose boolean indicating if timed logging is desired
#' @return something that's useful. TBD
#' @import veupathUtils
#' @import data.table
#' @export
rankedAbundance <- function(otu, method = c('median','max','q3','var'), cutoff=10, verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received OTU table with", NROW(otu), "samples and", (NCOL(otu)-1), "taxa."), verbose)

    # Reshape back to sample, taxonomicLevel, abundance
    formattedDT <- data.table::melt(otu, measure.vars=colnames(otu)[-1], variable.factor=F, variable.name='TaxonomicLevel', value.name="Abundance")

    rankedTaxa <- rankTaxa(formattedDT, method)

    # Extract top N taxa
    topN <- rankedTaxa[Abundance > 0, TaxonomicLevel]
    isCutoff <- FALSE
    if (NROW(topN) > cutoff) {
      topN <- topN[1:cutoff]
      computeMessage <- "Applied cutoff"
      isCutoff <- TRUE
    }

    keepCols <- c("SampleID", topN)
    dt = otu[, ..keepCols]

    veupathUtils::logWithTime("Finished ranking taxa", verbose)
    
    data.table::setnames(dt,'SampleID','record')
    

    # Collect attributes
    attr <- list('computationDetails' = computeMessage,
                 'parameters' = method,
                 'isCutoff' = isCutoff)
    
    #### Make into a function? Need to get entity from variables
    attr$computedVariableDetails <- list('id' = names(dt[, -c('record')]),
                                         'entity' = 'entity',
                                         'defaultRange' = c(0,1),
                                         'isCollection' = TRUE)

    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    return(dt)
}

#' Ranked abundance app
#'
#' This function returns the name of a json file with ranked abundance results.
#' 
#' @param otu data.frame with samples as rows, taxa as columns.
#' @param methods vector of strings indicating ranking methods to use. Must be a subset of c('median','max','q3','var').
#' @param cutoff integer indicating the maximium number of taxa to be kept after ranking.
#' @param verbose boolean indicating if timed logging is desired.
#' @return TBD
#' @import veupathUtils
#' @export
rankedAbundanceApp <- function(otu, methods=c('median','max','q3','var'), cutoff=10, verbose=c(TRUE, FALSE)) {

    verbose <- veupathUtils::matchArg(verbose)

    allMethods <- c('median','max','q3','var')
    # Allow any number of methods to be inputted
    if (is.null(methods)) methods <- allMethods
    if (!all(methods %in% allMethods)) {
      stop("Unaccepted method found in 'methods' argument. 'methods' must be a subset of c('median', 'max', 'q3','var').")
    }

    appResults <- lapply(methods, rankedAbundance, otu=otu, cutoff=cutoff, verbose=verbose)

    # Write to json file - debating whether to keep this in here or move elsewhere. Makes testing easier
    # outFileName <- writeAppResultsToJson(appResults, 'rankedAbundance', verbose = verbose)

    return(appResults)

}