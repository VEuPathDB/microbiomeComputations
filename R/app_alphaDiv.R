#' Alpha diversity
#'
#' This function returns alpha diversity values for each sample.
#' 
#' @param otu data.frame with samples as rows, taxa as columns
#' @param method string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'
#' @param verbose boolean indicating if timed logging is desired
#' @return something that's useful. TBD
#' @export
alphaDiv <- function(otu, method = c('shannon','simpson','evenness'), verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- plot.data::matchArg(method)
    verbose <- plot.data::matchArg(verbose)

    computeMessage <- ''
    plot.data::logWithTime(paste("Received OTU table with", NROW(otu), "samples and", (NCOL(otu)-1), "taxa."), verbose)
    

    # Compute alpha diversity
    if (identical(method, 'shannon') | identical(method, 'simpson')){

      alphaDivDT <- try(vegan::diversity(otu[, -c('SampleID')], method))
      computedVarLabel <- stringr::str_to_title(method)

    } else if (identical(method, 'evenness')) {

      alphaDivDT <- try(vegan::diversity(otu[, -c('SampleID')], 'shannon') / log(vegan::specnumber(otu)))
      computedVarLabel <- "Pielou\'s Evenness"
    }

    if (is.error(alphaDivDT)) {
      computeMessage <- paste('Error: alpha diversity', method, 'failed')
      # Also handle dt, results? Or just stop?
    } else {
      computeMessage <- paste('Computed', method, 'alpha diversity.')
    }

    # Assemble data table
    dt <- data.table('SampleID'= otu[['SampleID']],
                      'alphaDiv' = alphaDivDT)
    data.table::setnames(dt, 'alphaDiv', method)

    results <- list(
      'dt' = dt,
      'computedVariables'= names(dt[, -c('SampleID')]),
      'computedVariableLabels'= computedVarLabel,
      'computedAxisLabel' = jsonlite::unbox('Alpha Diversity'),
      'defaultRange' = c(0, 1),
      'computeDetails' = jsonlite::unbox(computeMessage),
      'computeName' = jsonlite::unbox(method)
    )

    plot.data::logWithTime(paste('Alpha diversity calculations completed with parameters method=', method, ', verbose =', verbose), verbose)
    return(results)
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
alphaDivApp <- function(otu, verbose = c(TRUE, FALSE)) {

    verbose <- plot.data::matchArg(verbose)

    methods <- c('shannon','simpson','evenness')

    computeResults <- lapply(methods, alphaDiv, otu=otu, verbose=verbose)


    # Merge all data into one data table and remove from computeResults
    dtList <- lapply(computeResults, function(x) return(x$dt))
    dt <- Reduce(function(...) merge(..., all = TRUE, by='SampleID'), dtList)
    computeResults <- lapply(computeResults, function(x) {x$dt <- NULL; return(x)})



    appResults <- list("data" = dt,
                        "stats" = NULL,
                        "metadata" = computeResults
    )

    # Write to json - nope, the following is in RServe general utils
    # outFileName <- writeListToJson(appResults, 'AlphaDiv')
    return(appResults)
}
