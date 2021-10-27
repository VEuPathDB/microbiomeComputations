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
    pcoa <- PCOA(dist)
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
    percentVar <- percentVar[1:k]

    # Perhaps this would be more clear if pcoaVar and computeDetails were attributes of dt? Worth making a class?
    results <- list(
      'dt' = dt,
      'computedVariables' = names(dt[, -c('SampleID')]),
      'computeName' = jsonlite::unbox(method),
      'computeDetails' = jsonlite::unbox(computeMessage),
      'pcoaVariance' = percentVar
    )

    plot.data::logWithTime(paste('Beta diversity calculations completed with parameters method =', method, ', k =', k, ', verbose =', verbose), verbose)
    return(results)
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

    computeResults <- betaDiv(otu, method, k, verbose)

    outDT <- computeResults$dt

    appResults <- list("data" = computeResults$dt,
                      "stats" = list(pcoaVariance = computeResults$pcoaVariance))

    computeResults$dt <- NULL
    computeResults$pcoaVariance <- NULL
    appResults$metadata <- computeResults


    # Write to json -- this lives in RServe utils so it's a no go here
    # outFileName <- writeListToJson(appResults, 'BetaDiv')
    return(appResults)

}