#' Alpha diversity
#'
#' This function returns alpha diversity values for each sample.
#' 
#' @param df data.frame with samples as rows, taxa as columns
#' @param sampleIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param method string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'
#' @param verbose boolean indicating if timed logging is desired
#' @return something that's useful. TBD
#' @importFrom vegan diversity
#' @importFrom stringi stri_trans_totitle
#' @import veupathUtils
#' @import data.table
#' @export
alphaDiv <- function(df, sampleIdColumn, method = c('shannon','simpson','evenness'), verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", NROW(df), "samples and", (NCOL(df)-1), "taxa."), verbose)
    

    # Compute alpha diversity
    if (identical(method, 'shannon') | identical(method, 'simpson')){

      alphaDivDT <- try(vegan::diversity(df[, -..sampleIdColumn], method))
      computedVarLabel <- stringi::stri_trans_totitle(method)

    } else if (identical(method, 'evenness')) {

      alphaDivDT <- try(vegan::diversity(df[, -..sampleIdColumn], 'shannon') / log(vegan::specnumber(df)))
      computedVarLabel <- "Pielou\'s Evenness"
    }

    if (veupathUtils::is.error(alphaDivDT)) {
      computeMessage <- paste('Error: alpha diversity', method, 'failed')
      # Also handle dt, results? Or just fail?
    } else {
      computeMessage <- paste('Computed', method, 'alpha diversity.')
    }

    # Assemble data table
    dt <- data.table('SampleID'= df[[sampleIdColumn]],
                      'alphaDiv' = alphaDivDT)

    data.table::setnames(dt, c(sampleIdColumn,'alphaDiversity'))
    
    # Collect attributes
    attr <- list('computationDetails' = computeMessage,
                 'parameters' = method,
                 'recordVariable' = sampleIdColumn)
    
    #### Make into a function? Need to get entity from variables
    attr$computedVariableDetails <- list('id' = names(dt[, -..sampleIdColumn]),
                                         'entity' = 'entity',
                                         'displayLabel' = computedVarLabel,
                                         'defaultRange' = c(0,1))

    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    veupathUtils::logWithTime(paste('Alpha diversity computation completed with parameters sampleIdColumn=', sampleIdColumn, ', method=', method, ', verbose =', verbose), verbose)
    
    return(dt)
}


#' Alpha diversity
#'
#' This function returns json with alpha diversity values for each sample.
#' 
#' @param df data.frame with samples as rows, taxa as columns
#' @param sampleIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param methods vector of strings defining the the beta diversity dissimilarity methods to use. Must be a subset of c('shannon','simpson','evenness').
#' @param verbose boolean indicating if timed logging is desired
#' @return something that's useful. TBD
#' @export
#' @import data.table
#' @import veupathUtils
alphaDivApp <- function(df, sampleIdColumn, methods = c('shannon','simpson','evenness'), verbose = c(TRUE, FALSE)) {

    verbose <- veupathUtils::matchArg(verbose)

    # Allow any number of methods to be inputted
    allMethods <- c('shannon','simpson','evenness')
    if (is.null(methods)) methods <- allMethods
    if (!all(methods %in% allMethods)) {
      stop("Unaccepted method found in 'methods' argument. 'methods' must be a subset of c('shannon','simpson','evenness').")
    }
    
    # Check that incoming df meets requirements
    if (!'data.table' %in% class(df)) {
      data.table::setDT(df)
    }
    if (!sampleIdColumn %in% names(df)) {
      stop("sampleIdColumn must exist as a column in df")
    }
    if (!all(unlist(lapply(df[, -..sampleIdColumn], is.numeric)))) {
      stop("All columns except the sampleIdColumn must be numeric")
    }

    appResults <- lapply(methods, alphaDiv, df=df, sampleIdColumn=sampleIdColumn, verbose=verbose)

    # Write to json file - debating whether to keep this in here or move elsewhere. Makes testing easier
    # outFileName <- writeAppResultsToJson(appResults, 'alphaDiv', verbose = verbose)

    return(appResults)
}
