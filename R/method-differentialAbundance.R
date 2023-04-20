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
  function(data, comparisonVariable, groupA, groupB, method = c('DESeq', 'ANCOMBC'), verbose = c(TRUE, FALSE)) standardGeneric("differentialAbundance"),
  signature = c("data")
)

#'@export
setMethod("differentialAbundance", signature("AbsoluteAbundanceData"), function(data, comparisonVariable, groupA, groupB, method = c('DESeq', 'ANCOMBC'), verbose = c(TRUE, FALSE)) {
    df <- data@data
    recordIdColumn <- data@recordIdColumn
    naToZero <- data@imputeZero
    ancestorIdColumns <- data@ancestorIdColumns
    allIdColumns <- c(recordIdColumn, ancestorIdColumns)
    sampleMetadata <- data@sampleMetadata


    # Initialize and check inputs
    method <- veupathUtils::matchArg(method)
    verbose <- veupathUtils::matchArg(verbose)

    # Check that incoming df meets requirements
    if (!'data.table' %in% class(df)) {
      data.table::setDT(df)
    }

    # Check that incoming metadata meets requirements
    if (!'data.table' %in% class(sampleMetadata)) {
      data.table::setDT(sampleMetadata)
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

    # Transpose abundance data to get a counts matrix with taxa on rows and samples as columns
    cts <- data.table::transpose(df[, -..recordIdColumn])
    rownames(cts) <- names(df[, -..recordIdColumn])
    colnames(cts) <- df[[recordIdColumn]]

    # Format metadata
    rownames(sampleMetadata) <- sampleMetadata[[recordIdColumn]]
    sampleMetadata <- sampleMetadata[, -..recordIdColumn]

    # Compute differential abundance
    if (identical(method, 'DESeq')) {


      # CHECK to ensure samples are in the same order for the cts and metadata dfs
      
      # Convert comparisonVariable column to factors so deseq doesn't have to do it

      # Create DESeqDataSet
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                            colData = sampleMetadata,
                                            design = as.formula(paste0("~",comparisonVariable)),
                                            tidy = FALSE) # consider changing to true so dont have to format metadata

      # Estimate size factors before running deseq to avoid errors about 0 counts
      geoMeans = apply(DESeq2::counts(dds), 1, function(x){exp(sum(log(x[x > 0]), na.rm=T) / length(x))})
      dds <- DESeq2::estimateSizeFactors(dds,geoMeans=geoMeans) # Alternatively let type ='iterate'. Pros/cons?

      # Run DESeq
      deseq_output <- DESeq2::DESeq(dds)

      # Extract results
      deseq_results <- DESeq2::results(deseq_output)

      # Format results for us
      statistics <- data.frame(log2foldChange = deseq_results$log2FoldChange,
                               pValue = deseq_results$pvalue,
                               adjustedPValue = deseq_results$padj,
                               pointID = rownames(cts))




    } else if (identical(method, 'ANCOMBC')) {

      se <- TreeSummarizedExperiment::TreeSummarizedExperiment(list(counts = cts), colData = sampleMetadata)

      # Currently getting this error: Error in is.infinite(o1) : default method not implemented for type 'list'
      # Ignoring for now.
      output_abs = ANCOMBC::ancombc2(data = se, assay_name = "counts", tax_level = NULL,
                  fix_formula = comparisonVariable, rand_formula = NULL,
                  p_adj_method = "holm", prv_cut=0,
                  group = comparisonVariable)

    } else {
      stop('Unaccepted differential abundance method. Accepted methods are DESeq and ANCOMBC.')
    }
    
    ## Make the result a Compute result. Return the samples that are left as the data. Return
    # the rest as statistics
    result <- new("ComputeResult")
    result@name <- 'differentialAbundance'
    result@recordIdColumn <- recordIdColumn
    result@ancestorIdColumns <- ancestorIdColumns
    result@statistics <- statistics


    # This should be the samples actually used in the method. Samples might get 
    # dropped because of whatever reasoning in the diff abund method.
    result@data <- df[, ..recordIdColumn]


    validObject(result) # not crazy about what i've done for ComputeResult. Considering
    # subclassing and making a ComputeResultWithStatistics subclass. Still have to navigate the computedVarMetadata though
    veupathUtils::logWithTime(paste('Differential abundance computation completed with parameters recordIdColumn=', recordIdColumn, ', method =', method, ', ..., verbose =', verbose), verbose)
    
    return(result)
})