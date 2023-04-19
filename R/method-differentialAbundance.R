#' Differential abundance
#'
#' This function returns the fold change and associated p value for a differential abundance analysis comparing samples in two groups.
#' 
#' @param data AbundanceData object
#' @param comparisonVariable string identifying the metadata column to be used when grouping samples.
#' @param groupA array of strings, each string matching a metadata value. Samples with this metadata value will be included in Group A
#' @param groupB array of strings, each string matching a metadata value. Samples with this metadata value will be included in Group B. Must be distinct from groupA.
#' @param method string defining the the differential abundance method. Accepted values are 'DESeq'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @useDynLib microbiomeComputations
#' @export
setGeneric("differentialAbundance",
  function(data, comparisonVariable, groupA, groupB, method = c('DESeq'), verbose = c(TRUE, FALSE)) standardGeneric("differentialAbundance"),
  signature = c("data")
)

#'@export
setMethod("differentialAbundance", signature("AbundanceData"), function(data, comparisonVariable, groupA, groupB, method = c('DESeq'), verbose = c(TRUE, FALSE)) {
    df <- data@data
    recordIdColumn <- data@recordIdColumn
    naToZero <- data@imputeZero
    ancestorIdColumns <- data@ancestorIdColumns
    allIdColumns <- c(recordIdColumn, ancestorIdColumns)
    metadata <- data@metadata


    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    # Check that incoming df meets requirements - consider moving to a validateOTU function or similar
    if (!'data.table' %in% class(df)) {
      data.table::setDT(df)
    }

    # Check that incoming metadata meets requirements - consider moving to a validateOTU function or similar
    if (!'data.table' %in% class(metadata)) {
      data.table::setDT(metadata)
    }

    computeMessage <- ''
    veupathUtils::logWithTime(paste("Received df table with", nrow(df), "samples and", (ncol(df)-1), "taxa."), verbose)

    if (naToZero) {
      # Replace NA values with 0
      veupathUtils::setNaToZero(df)
      veupathUtils::logWithTime("Replaced NAs with 0", verbose)
    }

    ### Take values in groupA, groupB, and turn them into a binary variable
    ## To do

    # Compute differential abundance
    if (identical(method, 'DESeq')) {

      # Lots of the following may get moved outside this conditional if it matches the prep for ancombc as well
      # Transpose abundance data to get a counts matrix with taxa on rows and samples as columns
      ### ANN set_rownames??
      cts <- data.table::transpose(df[, -..recordIdColumn])
      rownames(cts) <- names(df[, -..recordIdColumn])
      colnames(cts) <- df[[recordIdColumn]]


      # CHECK to ensure samples are in the same order for the cts and metadata dfs

      # Format metadata
      rownames(metadata) <- metadata[[recordIdColumn]]
      metadata <- metadata[, -..recordIdColumn]

      # Create DESeqDataSet
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                            colData = metadata,
                                            design = as.formula(paste("~",comparisonVariable)),
                                            tidy = FALSE) # consider changing to true so dont have to format metadata

      # Estimate size factors before running deseq to avoid errors about 0 counts
      geoMeans = apply(counts(dds), 1, function(x){exp(sum(log(x[x > 0]), na.rm=T) / length(x))})
      dds <- DESeq2::estimateSizeFactors(dds,geoMeans=geoMeans) # Alternatively let type ='iterate'. Pros/cons?

      # Run DESeq
      deseq_output <- DESeq2::DESeq(dds)

      # Extract results
      deseq_results <- DESeq2::results(deseq_output)

      # Format results for us
      statistics <- data.frame(log2FoldChange = deseq_results$log2FoldChange,
                               pvalue = deseq_results$pvalue,
                               adjustedPValue = deseq_results$padj,
                               pointID = rownames(cts))




    } else {
      stop('Unaccepted dissimilarity method. Accepted methods are bray, jaccard, and jsd.')
    }
    
    ## Make the result a Compute result. Return the samples that are left as the data. Return
    # the rest as statistics
    result <- new("ComputeResult")
    result@name <- 'differentialAbundance'
    result@recordIdColumn <- recordIdColumn
    result@ancestorIdColumns <- ancestorIdColumns
    result@otherSlot <- 'hi'
    result@statistics <- statistics

    # # Handle errors or return positive computeMessage
    # if (veupathUtils::is.error(dist)) {
    #   veupathUtils::logWithTime(paste('Differential abundance computation FAILED with parameters method=', method), verbose)
    #   stop() 
    # } else {
    #   veupathUtils::logWithTime("Computed dissimilarity matrix.", verbose)
    #   computeMessage <- paste(method, "dissimilarity matrix computation successful.")
    # }

    entity <- veupathUtils::strSplit(recordIdColumn,".", 4, 1)
    # result@computationDetails <- paste(computeMessage, ', pcoaVariance =', percentVar)
    # result@parameters <- paste('method =', method)   
    
    
    # makeVariableMetadataObject <- function(varName) {
    #   #bit hacky, see if you can think of something better
    #   plotRef <- ifelse(identical(varName, 'log2FoldChange'), 'xAxis', 'tooltip')

    #   veupathUtils::VariableMetadata(
    #              variableClass = veupathUtils::VariableClass(value = "computed"),
    #              variableSpec = veupathUtils::VariableSpec(variableId = varName, entityId = entity),
    #              plotReference = veupathUtils::PlotReference(value = plotRef),
    #             #  displayName = displayName, # Ann add logic for display name
    #              displayRangeMin = min(dt[[axisName]]),
    #              displayRangeMax = max(dt[[axisName]]),
    #              dataType = veupathUtils::DataType(value = "NUMBER"),
    #              dataShape = veupathUtils::DataShape(value = "CONTINUOUS")
    #   )
    # }

    ### This doesnot work yet...  
    # computedStatisticMetadata <- veupathUtils::VariableMetadataList(lapply(displayNames, makeVariableMetadataObject))

    # result@computedStatisticMetadata <- computedStatisticMetadata
    # names(dt) <- stripEntityIdFromColumnHeader(names(dt))
    result@data <- dt[[recordIdColumn]]

    validObject(result)
    veupathUtils::logWithTime(paste('Differential abundance computation completed with parameters recordIdColumn=', recordIdColumn, ', method =', method, ', ..., verbose =', verbose), verbose)
    
    return(result)
})