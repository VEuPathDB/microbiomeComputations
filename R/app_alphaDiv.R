#' Alpha diversity
#'
#' This function returns alpha diversity values for each sample.
#' 
#' @param df data.table with rows corresponding to records, one column indicating the record id, and all other columns corresponding to taxa abundance.
#' @param recordIdColumn string defining the name of the df column that specifies record ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param method string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'
#' @param verbose boolean indicating if timed logging is desired
#' @return data.table with columns recordIdColumn, "alphaDiversity".
#' @importFrom vegan diversity
#' @importFrom stringi stri_trans_totitle
#' @import veupathUtils
#' @import data.table
#' @export
alphaDiv <- function(df, recordIdColumn, method = c('shannon','simpson','evenness'), verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", NROW(df), "samples and", (NCOL(df)-1), "taxa."), verbose)
    

    # Compute alpha diversity
    if (identical(method, 'shannon') | identical(method, 'simpson')){

      alphaDivDT <- try(vegan::diversity(df[, -..recordIdColumn], method))
      computedVarLabel <- stringi::stri_trans_totitle(method)

    } else if (identical(method, 'evenness')) {

      alphaDivDT <- try(vegan::diversity(df[, -..recordIdColumn], 'shannon') / log(vegan::specnumber(df)))
      computedVarLabel <- "Pielou\'s Evenness"
    }

    if (veupathUtils::is.error(alphaDivDT)) {
      computeMessage <- paste('Error: alpha diversity', method, 'failed')
      # Also handle dt, results? Or just fail?
    } else {
      computeMessage <- paste('Computed', method, 'alpha diversity.')
    }

    # Assemble data table
    dt <- data.table('SampleID'= df[[recordIdColumn]],
                      'alphaDiv' = alphaDivDT)

    data.table::setnames(dt, c(recordIdColumn,'alphaDiversity'))
    
    # Collect attributes
    entity <- veupathUtils::strSplit(recordIdColumn,".", 4, 1)
    attr <- list('computationDetails' = computeMessage,
                 'parameters' = method,
                 'recordVariable' = recordIdColumn)
    
    #### Make into a function? Need to get entity from variables
    attr$computedVariableDetails <- list('id' = names(dt[, -..recordIdColumn]),
                                         'entity' = entity,
                                         'displayLabel' = computedVarLabel,
                                         'defaultRange' = c(0,1))
    # Add entity to column names
    data.table::setnames(dt, names(dt[, -..recordIdColumn]), paste0(entity,".",names(dt[, -..recordIdColumn])))

    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    veupathUtils::logWithTime(paste('Alpha diversity computation completed with parameters recordIdColumn=', recordIdColumn, ', method=', method, ', verbose =', verbose), verbose)
    
    return(dt)
}


#' Alpha diversity
#'
#' This function returns json with alpha diversity values for each sample.
#' 
#' @param df data.frame with samples as rows, taxa as columns
#' @param recordIdColumn string defining the name of the df column that specifies sample ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param methods vector of strings defining the the beta diversity dissimilarity methods to use. Must be a subset of c('shannon','simpson','evenness').
#' @param verbose boolean indicating if timed logging is desired
#' @return name of a jason file containing a list of data.tables, one for each method specified in methods. Each data.table contains columns recordIdColumn, "alphaDiversity", and an attribute "parameters" that records the method used.
#' @export
#' @import data.table
#' @import veupathUtils
alphaDivApp <- function(df, recordIdColumn, methods = c('shannon','simpson','evenness'), verbose = c(TRUE, FALSE)) {

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
    if (!recordIdColumn %in% names(df)) {
      stop("recordIdColumn must exist as a column in df")
    }
    if (!all(unlist(lapply(df[, -..recordIdColumn], is.numeric)))) {
      stop("All columns except the recordIdColumn must be numeric")
    }

    appResults <- lapply(methods, alphaDiv, df=df, recordIdColumn=recordIdColumn, verbose=verbose)

    # Write to json file - debating whether to keep this in here or move elsewhere. Makes testing easier
    # outFileName <- writeAppResultsToJson(appResults, 'alphaDiv', verbose = verbose)

    return(appResults)
}
