#' Differential abundance
#'
#' This function returns the fold change and associated p value for a differential abundance analysis comparing samples in two groups.
#' 
#' @param data AbsoluteAbundanceData object
#' @param comparisonVariable string identifying the metadata column to be used when grouping samples.
#' @param groupA array of strings, each string indicating a value in the comparisonVariable column. Samples with these metadata values will be included in Group A for the differential abundance calculation. Must have no overlap with groupB
#' @param groupB array of strings, each string indicating a value in the comparisonVariable column. Samples with these metadata values will be included in Group B for the differential abundance calculation. Must have no overlap with groupA.
#' @param method string defining the the differential abundance method. Accepted values are 'DESeq' and 'ANCOMBC'.
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @importFrom purrr none
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


    ## Initialize and check inputs
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

    ## Check that groups are provided, if needed, and if they are provided,
    ## that they match at least one value in the comparisonVariable column.
    print(sampleMetadata);
    print(names(sampleMetadata));
    uniqueComparisonVariableValues <- sort(unique(sampleMetadata[[comparisonVariable]]))
    
    if (!!length(groupA) && !!length(groupB)) {
      # Does each group contain at least one string that matches a value in the comparisonValue column?
      if (purrr::none(groupA,function(x, arr) { x %in% arr}, uniqueComparisonVariableValues)) {
        stop("At least one value in groupA must exist as a value in the comparisonValue sampleMetadata column.")
      }
      if (purrr::none(groupB,function(x, arr) { x %in% arr}, uniqueComparisonVariableValues)) {
        stop("At least one value in groupB must exist as a value in the comparisonValue sampleMetadata column.")
      }

      # Alert that we discard groupA/B values that are not found in the comparisonVariable
      if (any(!(groupA %in% uniqueComparisonVariableValues))) {
        veupathUtils::logWithTime("Found values in groupA that do not exist in the comparisonVariable. Removing these values.", verbose)
        groupA <- groupA[groupA %in% uniqueComparisonVariableValues]
      }
      if (any(!(groupB %in% uniqueComparisonVariableValues))) {
        veupathUtils::logWithTime("Found values in groupB that do not exist in the comparisonVariable. Removing these values.", verbose)
        groupB <- groupB[groupB %in% uniqueComparisonVariableValues]
      }

      # Do not allow duplicated values
      if (!!length(intersect(groupA, groupB))) {
        veupathUtils::logWithTime("groupA and groupB cannot share members.", verbose)
        stop()
      }

    } else if (length(uniqueComparisonVariableValues) == 2) {
      # Ignore any groups passed to us and set the groups to be these two values
      veupathUtils::logWithTime("Only two values found in the comparisonVariable sampleMetadata column. Using these two values as groupA and groupB.", verbose)
      groupA <- c(uniqueComparisonVariableValues[1])
      groupB <- c(uniqueComparisonVariableValues[2])
    } else {
      stop("Must supply two groups (groupA and groupB) for the differential abundance calculation.")
    }

    # Subset to only include samples with metadata defined in groupA and groupB
    sampleMetadata <- sampleMetadata[get(comparisonVariable) %in% c(groupA, groupB), ]
    keepSamples <- sampleMetadata[[recordIdColumn]]
    veupathUtils::logWithTime(paste0("Found ",length(keepSamples)," samples with ", comparisonVariable, "in either groupA or groupB. The calculation will continue with only these samples."), verbose)

    # Subset the original data based on the kept samples
    df <- df[get(recordIdColumn) %in% keepSamples, ]
    sampleMetadata <- sampleMetadata[get(recordIdColumn) %in% keepSamples, ]

    # Turn comparisonVariable into a binary variable
    sampleMetadata[get(comparisonVariable) %in% groupA, c(comparisonVariable)] <- 'groupA'
    sampleMetadata[get(comparisonVariable) %in% groupB, c(comparisonVariable)] <- 'groupB'


    ## Format data for the different differential abundance methods.

    # First, transpose abundance data to get a counts matrix with taxa as rows and samples as columns
    counts <- data.table::transpose(df[, -..recordIdColumn])
    rownames(counts) <- names(df[, -..recordIdColumn])
    colnames(counts) <- df[[recordIdColumn]]

    # Next, format metadata. Recall samples are rows and variables are columns
    rownames(sampleMetadata) <- sampleMetadata[[recordIdColumn]]

    # Finally, check to ensure samples are in the same order in counts and metadata. Both DESeq
    # and ANCOMBC expect the order to match, and will not perform this check.
    if (!identical(rownames(sampleMetadata), colnames(counts))){
      # Reorder sampleMetadata to match counts
      veupathUtils::logWithTime("Sample order differs between data and metadata. Reordering data based on the metadata sample order.", verbose)
      data.table::setcolorder(counts, rownames(sampleMetadata))
    }
    veupathUtils::logWithTime(paste0("Abundance data formatted for differential abundance computation. Proceeding with method=",method), verbose)
    
    ## Compute differential abundance
    if (identical(method, 'DESeq')) {

      deseq_output <- try({

        # Create DESeqDataSet (dds)
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                              colData = sampleMetadata,
                                              design = as.formula(paste0("~",comparisonVariable)),
                                              tidy = FALSE)

        # Estimate size factors before running deseq to avoid errors about 0 counts
        geoMeans = apply(DESeq2::counts(dds), 1, function(x){exp(sum(log(x[x > 0]), na.rm=T) / length(x))})
        dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geoMeans)

        # Run DESeq
        deseq_output <- DESeq2::DESeq(dds)
      })

      if (veupathUtils::is.error(deseq_output)) {
        veupathUtils::logWithTime(paste0('Differential abundance FAILED with parameters recordIdColumn=', recordIdColumn, ', method =', method, ', verbose =', verbose), verbose)
        stop()
      }

      # Extract deseq results
      deseq_results <- DESeq2::results(deseq_output)

      # Format results for easier access
      statistics <- data.frame(log2foldChange = deseq_results$log2FoldChange,
                               pValue = deseq_results$pvalue,
                               adjustedPValue = deseq_results$padj,
                               pointID = rownames(counts))


    } else if (identical(method, 'ANCOMBC')) {

      se <- TreeSummarizedExperiment::TreeSummarizedExperiment(list(counts = counts), colData = sampleMetadata)

      # Currently getting this error: Error in is.infinite(o1) : default method not implemented for type 'list'
      # Ignoring for now.
      output_abs = ANCOMBC::ancombc2(data = se, assay_name = "counts", tax_level = NULL,
                  fix_formula = comparisonVariable, rand_formula = NULL,
                  p_adj_method = "holm", prv_cut=0,
                  group = comparisonVariable)

    } else {
      stop('Unaccepted differential abundance method. Accepted methods are "DESeq" and "ANCOMBC".')
    }
    veupathUtils::logWithTime(paste0('Completed method=',method,'. Formatting results.'), verbose)
    
    ## Construct the ComputeResult
    result <- new("ComputeResult")
    result@name <- 'differentialAbundance'
    result@recordIdColumn <- recordIdColumn
    result@ancestorIdColumns <- ancestorIdColumns
    result@statistics <- statistics
    result@parameters <- paste0("comparisonVariable = ", comparisonVariable, ", groupA = ", groupA, ", groupB = ", groupB, ', method = ', method)


    # The resulting data should contain only the samples actually used.
    result@data <- df[, ..allIdColumns]
    names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))


    validObject(result)
    veupathUtils::logWithTime(paste('Differential abundance computation completed with parameters recordIdColumn=', recordIdColumn, "comparisonVariable = ", comparisonVariable, ", groupA = ", groupA, ", groupB = ", groupB, ', method = ', method), verbose)
    
    return(result)
})