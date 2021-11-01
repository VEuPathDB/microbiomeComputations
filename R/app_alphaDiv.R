#' Alpha diversity
#'
#' This function returns alpha diversity values for each sample.
#' 
#' @param otu data.frame with samples as rows, taxa as columns
#' @param method string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'
#' @param verbose boolean indicating if timed logging is desired
#' @return something that's useful. TBD
#' @importFrom vegan diversity
#' @importFrom stringi stri_trans_totitle
#' @import veupathUtils
#' @import data.table
#' @export
alphaDiv <- function(otu, method = c('shannon','simpson','evenness'), verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received OTU table with", NROW(otu), "samples and", (NCOL(otu)-1), "taxa."), verbose)
    

    # Compute alpha diversity
    if (identical(method, 'shannon') | identical(method, 'simpson')){

      alphaDivDT <- try(vegan::diversity(otu[, -c('SampleID')], method))
      computedVarLabel <- stringi::stri_trans_totitle(method)

    } else if (identical(method, 'evenness')) {

      alphaDivDT <- try(vegan::diversity(otu[, -c('SampleID')], 'shannon') / log(vegan::specnumber(otu)))
      computedVarLabel <- "Pielou\'s Evenness"
    }

    if (veupathUtils::is.error(alphaDivDT)) {
      computeMessage <- paste('Error: alpha diversity', method, 'failed')
      # Also handle dt, results? Or just fail?
    } else {
      computeMessage <- paste('Computed', method, 'alpha diversity.')
    }

    # Assemble data table
    dt <- data.table('SampleID'= otu[['SampleID']],
                      'alphaDiv' = alphaDivDT)
    data.table::setnames(dt, 'alphaDiv', method)

    attr <- list(
      'computedVariables'= names(dt[, -c('SampleID')]),
      'computedVariableLabels'= computedVarLabel,
      'computedAxisLabel' = 'Alpha Diversity',
      'defaultRange' = c(0, 1),
      'computeDetails' = computeMessage
    )
    
    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    veupathUtils::logWithTime(paste('Alpha diversity calculations completed with parameters method=', method, ', verbose =', verbose), verbose)
    
    return(dt)
}


#' Alpha diversity
#'
#' This function returns json with alpha diversity values for each sample.
#' 
#' @param otu data.frame with samples as rows, taxa as columns
#' @param verbose boolean indicating if timed logging is desired
#' @return something that's useful. TBD
#' @export
#' @import data.table
#' @import veupathUtils
alphaDivApp <- function(otu, verbose = c(TRUE, FALSE)) {

    verbose <- veupathUtils::matchArg(verbose)

    methods <- c('shannon','simpson','evenness')

    appResults <- lapply(methods, alphaDiv, otu=otu, verbose=verbose)
    
    names(appResults) <- methods

    # Write to json file - debating whether to keep this in here or move elsewhere. Makes testing easier
    # outFileName <- writeAppResultsToJson(appResults, 'alphaDiv', verbose = T)

    return(appResults)
}
