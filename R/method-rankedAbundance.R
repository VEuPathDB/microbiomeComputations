#' Ranked abundance
#'
#' This function returns abundances, ranked by a selected ranking function
#' 
#' @param data AbundanceData object
#' @param method string defining the ranking strategy by which to order the taxa. Accepted values are 'median','max','q3',and 'variance'. Note that taxa that return a value of 0 for a given method will not be included in the results.
#' @param cutoff integer indicating the maximium number of taxa to be kept after ranking.
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @importFrom S4Vectors SimpleList
#' @export
setGeneric("rankedAbundance",
  function(data, method = c('median','max','q3','variance'), cutoff=10, verbose = c(TRUE, FALSE)) standardGeneric("rankedAbundance"),
  signature = c("data")
)

#'@export 
setMethod("rankedAbundance", signature("AbundanceData"), function(data, method = c('median','max','q3','variance'), cutoff=10, verbose = c(TRUE, FALSE)) {
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

    # Reshape back to sample, taxonomicLevel, abundance
    formattedDT <- data.table::melt(df, measure.vars=colnames(df[, -..recordIdColumn]), variable.factor=F, variable.name='TaxonomicLevel', value.name="Abundance")

    rankedTaxa <- rankTaxa(formattedDT, method)

    # Extract top N taxa
    topN <- rankedTaxa[Abundance > 0, TaxonomicLevel]
    isCutoff <- FALSE
    if (length(topN) > cutoff) {
      topN <- topN[1:cutoff]
      computeMessage <- paste("Only returning top", cutoff, "results.")
      isCutoff <- TRUE
    }

    keepCols <- c(recordIdColumn, topN)
    dt = df[, ..keepCols]

    veupathUtils::logWithTime("Finished ranking taxa", verbose)
    
    result <- new("ComputeResult")
    result@name <- 'rankedAbundance'
    result@recordIdColumn <- recordIdColumn
    result@data <- dt

    entity <- veupathUtils::strSplit(recordIdColumn,".", 4, 1)
    result@computationDetails <- computeMessage
    # TODO this could be a named vector or something maybe instead
    # also, not sure isCutoff is a param, maybe put it in computationDetails?
    result@parameters <- paste0('method = ',method, ', isCutoff = ', isCutoff)

    collectionMemberVariableIds <- unlist(lapply(names(dt[, -..recordIdColumn]), veupathUtils::strSplit, ".", 4, 2))

    makeVariableSpecs <- function(variableId) {
	    veupathUtils::VariableSpec(variableId = variableId, entityId = entity)
    }

    computedVariableMetadata <- veupathUtils::VariableMetadata(
                 variableClass = veupathUtils::VariableClass(value = "computed"),
                 variableSpec = veupathUtils::VariableSpec(variableId = "rankedAbundance", entityId = entity),
                 plotReference = veupathUtils::PlotReference(value = "xAxis"),
                 displayName = "To be determined by client",
                 displayRangeMin = 0,
                 displayRangeMax = 1,
                 dataType = veupathUtils::DataType(value = "NUMBER"),
                 dataShape = veupathUtils::DataShape(value = "CONTINUOUS"),
                 isCollection = TRUE,
                 members = veupathUtils::VariableSpecList(S4Vectors::SimpleList(lapply(collectionMemberVariableIds, makeVariableSpecs)))
      )
    
    result@computedVariableMetadata <- veupathUtils::VariableMetadataList(S4Vectors::SimpleList(computedVariableMetadata))
    
    validObject(result)
    veupathUtils::logWithTime(paste('Ranked abundance computation completed with parameters recordIdColumn=', recordIdColumn, ', method =', method, ', cutoff =', cutoff, ', naToZero = ', naToZero, ', verbose =', verbose), verbose)

    return(result)
  })
