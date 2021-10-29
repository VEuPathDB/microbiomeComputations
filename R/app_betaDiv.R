#' Beta diversity
#'
#' This function returns pcoa coordinates calculated from the beta diversity dissimilarity matrix.
#' 
#' @param otu data.frame with samples as rows, taxa as columns
#' @param method string defining the the beta diversity dissimilarity method. Accepted values are 'bray','jaccard', and 'jsd'
#' @param k integer determining the number of pcoa dimensions to return
#' @param verbose boolean indicating if timed logging is desired
#' @return who knows
#' @export
betaDiv <- function(otu,
                    method = c('bray','jaccard','jsd'),
                    k = 2,
                    verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- plot.data::matchArg(method)
    verbose <- plot.data::matchArg(verbose)

    computeMessage <- ''
    plot.data::logWithTime(paste("Received OTU table with", NROW(otu), "samples and", (NCOL(otu)-1), "taxa."), verbose)

    # Compute beta diversity using given dissimilarity method
    if (identical(method, 'bray') | identical(method, 'jaccard')) {

      dist <- vegan::vegdist(otu[, -c('SampleID')], method=method, binary=TRUE)

    } else if (identical(method, 'jsd')) {

      otuMat <- matrix(as.numeric(unlist(otu[, -c("SampleID")])), nrow=NROW(otu))
      dist <- jsd(t(otuMat))
      dist <- as.dist(dist)

    } else {
      stop('Unaccepted dissimilarity method. Accepted methods are bray, jaccard, and jsd.')
    }
    plot.data::logWithTime("Computed dissimilarity matrix.", verbose)

    # Ordination
    ## Need to handle how this might err
    pcoa <- ape::pcoa(dist)
    dt <- data.table::as.data.table(pcoa$vectors)
    computeMessage <- paste("PCoA returned results for", NCOL(dt), "dimensions.")

    dt$SampleID <- otu[['SampleID']]
    data.table::setcolorder(dt, c('SampleID'))
    plot.data::logWithTime("Finished ordination step.", verbose)

    # Extract percent variance
    eigenvecs <- pcoa$values$Relative_eig
    percentVar <- round(100*(eigenvecs / sum(eigenvecs)), 1)

    # Keep dims 1:k
    # We should keep the same number of percentVar values as cols in the data table. However, i think we're letting the user download lots of columns? So perhaps we shouldn't have k at all? A plot can use however many it needs.
    # For now returning data and percentVar for how much is in the plot.
    percentVar <- percentVar[1:k]
    keepCols <- c('SampleID',names(dt)[2:(k+1)])
    dt <- dt[, ..keepCols]

    #### Need to add back computed Variable Labels
    attr <- list(
      'computedVariables' = names(dt[, -c('SampleID')]),
      'computeDetails' = computeMessage,
      'pcoaVariance' = percentVar
    )
    
    plot.data::setAttrFromList(dt, attr, removeExtraAttrs = F)

    plot.data::logWithTime(paste('Beta diversity calculations completed with parameters method =', method, ', k =', k, ', verbose =', verbose), verbose)
    
    return(dt)
}

#' Beta diversity app
#'
#' This function returns the name of a json file with beta diversity results.
#' 
#' @param otu data.frame with samples as rows, taxa as columns
#' @param method string defining the the beta diversity dissimilarity method. Accepted values are 'bray','jaccard', and 'jsd'
#' @param k integer determining the number of pcoa dimensions to return
#' @param verbose boolean indicating if timed logging is desired
#' @return we'll see.
#' @export
betaDivApp <- function(otu,
                    method = c('bray','jaccard','jsd'),
                    k = 2,
                    verbose = c(TRUE, FALSE)) {

    otu <- data.table::setDT(otu)

    method <- plot.data::matchArg(method)
    verbose <- plot.data::matchArg(verbose)

    appResults <- lapply(method, betaDiv, otu=otu, k=k, verbose=verbose)
    
    names(appResults) <- method

    # Write to json file - debating whether to keep this in here or move elsewhere. Makes testing easier
    # outFileName <- writeAppResultsToJson(appResults, 'betaDiv', verbose = T)

    return(appResults)

}