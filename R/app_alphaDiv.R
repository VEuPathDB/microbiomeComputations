#' Alpha diversity
#'
#' This function returns alpha diversity values for each sample.
#' 
#' @param df data.table with rows corresponding to records, one column indicating the record id, and all other columns corresponding to taxa abundance.
#' @param recordIdColumn string defining the name of the df column that specifies record ids. Note, all other columns must be numeric and will be treated as abundance values.
#' @param method string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'
#' @param naToZero boolean indicating if NAs in numeric columns should be replaced with 0
#' @param verbose boolean indicating if timed logging is desired
#' @return data.table with columns recordIdColumn, "alphaDiversity".
#' @importFrom vegan diversity
#' @importFrom stringi stri_trans_totitle
#' @import veupathUtils
#' @import data.table
#' @export
alphaDiv <- function(df, recordIdColumn, method = c('shannon','simpson','evenness'), naToZero = c(TRUE, FALSE), verbose = c(TRUE, FALSE)) {

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    naToZero <- veupathUtils::matchArg(naToZero)
    verbose <- veupathUtils::matchArg(verbose)

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
    
    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", nrow(df), "samples and", (ncol(df)-1), "taxa."), verbose)
    
    if (naToZero) {
      # Replace NA values with 0
      veupathUtils::setNaToZero(df)
      veupathUtils::logWithTime("Replaced NAs with 0", verbose)
    }

    # Compute alpha diversity
    if (identical(method, 'shannon') | identical(method, 'simpson')){

      alphaDivDT <- try(vegan::diversity(df[, -..recordIdColumn], method))
      computedVarLabel <- paste(stringi::stri_trans_totitle(method), 'Diversity')

    } else if (identical(method, 'evenness')) {

      alphaDivDT <- try(vegan::diversity(df[, -..recordIdColumn], 'shannon') / log(vegan::specnumber(df)))
      computedVarLabel <- "Pielou\'s Evenness"
    }

    # Handle errors or return positive computeMessage
    if (veupathUtils::is.error(alphaDivDT)) {
      
      computeMessage <- paste('Error: alpha diversity', method, 'failed:', attr(alphaDivDT,'condition')$message)
      
      # Return only recordIdColumn and expected attributes
      dt <- df[, ..recordIdColumn]
      
      attr <- list('computationDetails' = computeMessage,
                   'parameters' = character())
      
      computedVariableDetails <- list('variableId' = character(),
                                       'entityId' = character(),
                                       'dataType' = character(),
                                       'dataShape' = character(),
                                       'values' = numeric())
      
      computedVariableMetadata <- list('displayName' = character(),
                                       'displayRangeMin' = character(),
                                       'displayRangeMax' = character())
      
      attr$computedVariable <- list('computedVariableDetails' = computedVariableDetails,
                                     'computedVariableMetadata' = computedVariableMetadata)
      
      veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)
      veupathUtils::logWithTime(paste('Alpha diversity computation FAILED with parameters recordIdColumn=', recordIdColumn, ', method=', method, ', naToZero = ', naToZero, ', verbose =', verbose), verbose)
      
      return(dt)
      
    } else {
      computeMessage <- paste('Computed', method, 'alpha diversity.')
      veupathUtils::logWithTime(paste(method, 'alpha diversity computation complete.'), verbose)
    }

    # Assemble data table
    dt <- data.table('SampleID'= df[[recordIdColumn]],
                      'alphaDiv' = alphaDivDT)

    data.table::setnames(dt, c(recordIdColumn,'alphaDiversity'))
    
    # Collect attributes
    entity <- veupathUtils::strSplit(recordIdColumn, ".", 4, 1)
    attr <- list('computationDetails' = computeMessage,
                 'parameters' = method)

    computedVariableMetadata <- new("VariableMetadata",
                 variableClass = new("VariableClass", value = "computed"),
                 variableSpec = new("VariableSpec", variableId = names(dt[, -..recordIdColumn]), entityId = entity),
                 displayName = computedVarLabel,
                 displayRangeMin = 0,
                 displayRangeMax = 1,
                 dataType = new("DataType", value = "NUMBER"),
                 dataShape = new("DataShape", value = "CONTINUOUS")
      )
      
    attr$computedVariable <- computedVariableMetadata
    
    # Add entity to column names
    data.table::setnames(dt, names(dt[, -..recordIdColumn]), paste0(entity,".",names(dt[, -..recordIdColumn])))

    veupathUtils::setAttrFromList(dt, attr, removeExtraAttrs = F)

    veupathUtils::logWithTime(paste('Alpha diversity computation completed with parameters recordIdColumn=', recordIdColumn, ', method=', method, ', verbose =', verbose), verbose)
    
    return(dt)
}
