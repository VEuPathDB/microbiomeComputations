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

    # Compute differential abundance
    if (identical(method, 'DESeq')) {

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
                                            tidy = FALSE)

      # Estimate size factors before running deseq to avoid errors about 0 counts
      # NEEDS TO MOVE calculate geometric means prior to estimate size factors
      gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }
      geoMeans = apply(counts(dds), 1, gm_mean)
      dds <- DESeq2::estimateSizeFactors(dds,geoMeans=geoMeans)

      # Run DESeq
      deseq_output <- DESeq2::DESeq(dds)

      # Extract results
      deseq_results <- DESeq2::results(deseq_output)



    } else {
      stop('Unaccepted dissimilarity method. Accepted methods are bray, jaccard, and jsd.')
    }
    
    result <- new("ComputeResult")
    result@name <- 'differentialAbundance'
    # result@recordIdColumn <- recordIdColumn May be pointID now?
    # result@ancestorIdColumns <- ancestorIdColumns  needed?

    # # Handle errors or return positive computeMessage
    # if (veupathUtils::is.error(dist)) {
    #   veupathUtils::logWithTime(paste('Differential abundance computation FAILED with parameters method=', method), verbose)
    #   stop() 
    # } else {
    #   veupathUtils::logWithTime("Computed dissimilarity matrix.", verbose)
    #   computeMessage <- paste(method, "dissimilarity matrix computation successful.")
    # }

    # Ordination
    pcoa <- ape::pcoa(dist)
    dt <- data.table::as.data.table(pcoa$vectors)
    dt <- dt[, 1:k]
    # Remove dots from names
    data.table::setnames(dt, stringi::stri_replace_all_fixed(names(dt),".",""))
    computeMessage <- paste(computeMessage, "PCoA returned results for", ncol(dt), "dimensions.")

    dt <- cbind(dt, df[, ..allIdColumns])
    data.table::setcolorder(dt, allIdColumns)
    veupathUtils::logWithTime("Finished ordination step.", verbose)

    # Extract percent variance
    eigenvecs <- pcoa$values$Relative_eig
    percentVar <- round(100*(eigenvecs / sum(eigenvecs)), 1)

    # Keep dims 1:k
    # We should keep the same number of percentVar values as cols in the data table. However, i think we're letting the user download lots of columns? So perhaps we shouldn't have k at all? A plot can use however many it needs.
    # For now returning data and percentVar for how much is in the plot.
    percentVar <- percentVar[1:k]
    
    entity <- veupathUtils::strSplit(recordIdColumn,".", 4, 1)
    result@computationDetails <- paste(computeMessage, ', pcoaVariance =', percentVar)
    result@parameters <- paste('method =', method)   
    
    axesNames <- names(dt[, -..allIdColumns])
    displayNames <- paste0(axesNames, " ", sprintf(percentVar,fmt = '%#.1f'), "%")

    makeVariableMetadataObject <- function(displayName) {
      axisName <- veupathUtils::strSplit(displayName, " ")
      #bit hacky, see if you can think of something better
      plotRef <- ifelse(grepl('Axis1', displayName, fixed=T), 'xAxis', 'yAxis')

      veupathUtils::VariableMetadata(
                 variableClass = veupathUtils::VariableClass(value = "computed"),
                 variableSpec = veupathUtils::VariableSpec(variableId = axisName, entityId = entity),
                 plotReference = veupathUtils::PlotReference(value = plotRef),
                 displayName = displayName,
                 displayRangeMin = min(dt[[axisName]]),
                 displayRangeMax = max(dt[[axisName]]),
                 dataType = veupathUtils::DataType(value = "NUMBER"),
                 dataShape = veupathUtils::DataShape(value = "CONTINUOUS")
      )
    }
          
    computedVariableMetadata <- veupathUtils::VariableMetadataList(lapply(displayNames, makeVariableMetadataObject))

    result@computedVariableMetadata <- computedVariableMetadata
    names(dt) <- stripEntityIdFromColumnHeader(names(dt))
    result@data <- dt

    validObject(result)
    veupathUtils::logWithTime(paste('Beta diversity computation completed with parameters recordIdColumn=', recordIdColumn, ', method =', method, ', k =', k, ', verbose =', verbose), verbose)
    
    return(result)
})