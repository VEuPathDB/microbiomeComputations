#' Alpha diversity
#'
#' This function returns alpha diversity values for each sample.
#' 
#' @param data AbundanceData object
#' @param method string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputedResult object
#' @importFrom vegan diversity
#' @importFrom stringi stri_trans_totitle
#' @import veupathUtils
#' @import data.table
#' @export
setGeneric("alphaDiv",
  function(data, method, verbose) standardGeneric("alphaDiv"),
  signature = c("data")
)

#'@export 
setMethod("alphaDiv", signature("AbundanceData"), function(data, method = c('shannon','simpson','evenness'), verbose = c(TRUE, FALSE)) {
    df <- data@data
    recordIdColumn <- data@recordIdColumn
    naToZero <- data@imputeZero

    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    # Check that incoming df meets requirements
    if (!'data.table' %in% class(df)) {
      data.table::setDT(df)
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

    result <- new("ComputedResult")

    # Handle errors or return positive computeMessage
    if (veupathUtils::is.error(alphaDivDT)) {
      
      computeMessage <- paste('Error: alpha diversity', method, 'failed:', attr(alphaDivDT,'condition')$message)
      
      # Return only recordIdColumn and expected attributes
      dt <- df[, ..recordIdColumn]
      
      result@computationDetails = computeMessage
      result@computedVariableMetadata <- veupathUtils::VariableMetadataList(S4Vectors::SimpleList(veupathUtils::VariableMetadata()))
      
      veupathUtils::logWithTime(paste('Alpha diversity computation FAILED with parameters, method=', method), verbose)
      
      return(dt)
      
    } else {
      computeMessage <- paste('Computed', method, 'alpha diversity.')
      veupathUtils::logWithTime(paste(method, 'alpha diversity computation complete.'), verbose)
    }

    # Assemble data table
    dt <- data.table('SampleID'= df[[recordIdColumn]],
                      'alphaDiv' = alphaDivDT)

    data.table::setnames(dt, c(recordIdColumn,'alphaDiversity'))

    entity <- veupathUtils::strSplit(recordIdColumn, ".", 4, 1)
    result@computationDetails <- computeMessage
    result@parameters <- method

    computedVariableMetadata <- veupathUtils::VariableMetadata(
                 variableClass = veupathUtils::VariableClass(value = "computed"),
                 variableSpec = veupathUtils::VariableSpec(variableId = names(dt[, -..recordIdColumn]), entityId = entity),
                 plotReference = veupathUtils::PlotReference(value = "yAxis"),
                 displayName = computedVarLabel,
                 displayRangeMin = 0,
                 displayRangeMax = 1,
                 dataType = veupathUtils::DataType(value = "NUMBER"),
                 dataShape = veupathUtils::DataShape(value = "CONTINUOUS")
      )
      
    result@computedVariableMetadata <- veupathUtils::VariableMetadataList(S4Vectors::SimpleList(computedVariableMetadata))
    
    # Add entity to column names
    data.table::setnames(dt, names(dt[, -..recordIdColumn]), paste0(entity,".",names(dt[, -..recordIdColumn])))
    result@data <- dt 
    
    validObject(result)
    veupathUtils::logWithTime(paste('Alpha diversity computation completed with parameters method=', method), verbose)
    
    return(result)
})
