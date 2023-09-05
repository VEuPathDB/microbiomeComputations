#' Differential abundance
#'
#' This function returns the fold change and associated p value for a differential abundance analysis comparing samples in two groups.
#' 
#' @param data AbsoluteAbundanceData object
#' @param comparator Comparator object specifying the variable and values or bins to be used in dividing samples into groups.
#' @param method string defining the the differential abundance method. Accepted values are 'DESeq' and 'ANCOMBC'.
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @import veupathUtils
#' @import data.table
#' @importFrom purrr none
#' @importFrom purrr discard
#' @useDynLib microbiomeComputations
#' @export
setGeneric("differentialAbundance",
  function(data, comparator, method = c('DESeq', 'ANCOMBC'), verbose = c(TRUE, FALSE)) standardGeneric("differentialAbundance"),
  signature = c("data", "comparator")
)

#'@export
setMethod("differentialAbundance", signature("AbsoluteAbundanceData", "Comparator"), function(data, comparator, method = c('DESeq', 'ANCOMBC'), verbose = c(TRUE, FALSE)) {
    df <- data@data
    recordIdColumn <- data@recordIdColumn
    naToZero <- data@imputeZero
    ancestorIdColumns <- data@ancestorIdColumns
    allIdColumns <- c(recordIdColumn, ancestorIdColumns)
    sampleMetadata <- copy(data@sampleMetadata)
    comparatorColName <- veupathUtils::getColName(comparator@variable@variableSpec)


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

    # Replace NA values with 0
    if (naToZero) {
      veupathUtils::setNaToZero(df)
      veupathUtils::logWithTime("Replaced NAs with 0", verbose)
    }


    # Subset to only include samples with metadata defined in groupA and groupB
    if (identical(comparator@variable@dataShape@value, "CONTINUOUS")) {
      # Ensure bin starts and ends are numeric
      comparator@groupA <- veupathUtils::BinList(lapply(comparator@groupA, function(bin) {
        return(veupathUtils::Bin(
          binStart = as.numeric(bin@binStart),
          binEnd = as.numeric(bin@binEnd),
          binLabel = bin@binLabel
        ))
      }))
      comparator@groupB <- veupathUtils::BinList(lapply(comparator@groupB, function(bin) {
        return(veupathUtils::Bin(
          binStart = as.numeric(bin@binStart),
          binEnd = as.numeric(bin@binEnd),
          binLabel = bin@binLabel
        ))
      }))
      # We need to turn the numeric comparison variable into a categorical one with those values
      # that fall within group A or group B bins marked with some string we know.

      # Collect all instances where the comparatorColName has values in the bins from each group.
      # So inGroupA is a vector with 0 if the value in comparatorColName is not within any of the group A bins and >0 otherwise.
      lapply(comparator@groupA, print)
      lapply(comparator@groupB, print)
      print(sampleMetadata[[comparatorColName]])
      inGroupA <- veupathUtils::whichValuesInBinList(sampleMetadata[[comparatorColName]], comparator@groupA)
      inGroupB <- veupathUtils::whichValuesInBinList(sampleMetadata[[comparatorColName]], comparator@groupB)

      print(inGroupA)
      print(inGroupB)

      # Eventually move this check to Comparator validation. See #47
      if ((any(inGroupA * inGroupB) > 0)) {
        stop("Group A and Group B cannot have overlapping bins.")
      }

      # Make the comparatorColName a character vector and replace the in-group values with a bin.
      sampleMetadata[, (comparatorColName) := as.character(get(comparatorColName))]
      
      # Now we can reassign groupA and groupB and can replace values in sampleMetadata our new group values
      # We don't care about the values in the comparisonVariable column anymore. They were only
      # useful to help us assign groups.
      sampleMetadata[inGroupA, c(comparatorColName)] <- "groupA"
      sampleMetadata[inGroupB, c(comparatorColName)] <- "groupB"

      # Finally, subset the sampleMetadata to only include those samples in groupA or B
      sampleMetadata <- sampleMetadata[get(comparatorColName) %in% c("groupA", "groupB"), ]

    } else {
      # The comparator must be ordinal, binary, or categorical
      groupAValues <- getGroupLabels(comparator, "groupA")
      groupBValues <- getGroupLabels(comparator, "groupB")

      # Filter sampleMetadata to keep only those samples that are labeled as groupA or groupB. Filter
      # data *before* reassigning values to 'groupA' and 'groupB' to avoid issues with the original variable
      # value being 'groupA' or 'groupB'
      sampleMetadata <- sampleMetadata[get(comparatorColName) %in% c(groupAValues, groupBValues), ]

      # Turn comparatorColName into a binary variable
      sampleMetadata[get(comparatorColName) %in% groupAValues, c(comparatorColName)] <- 'groupA'
      sampleMetadata[get(comparatorColName) %in% groupBValues, c(comparatorColName)] <- 'groupB'
    } 

    # sampleMetadata has already been filtered so it now only contains the samples we care about
    keepSamples <- sampleMetadata[[recordIdColumn]]
    if (!length(keepSamples)) {
      stop("No samples remain after subsetting based on the comparator variable.")
    }
    veupathUtils::logWithTime(paste0("Found ",length(keepSamples)," samples with ", comparatorColName, "in either groupA or groupB. The calculation will continue with only these samples."), verbose)

    # Subset the abundance data based on the kept samples
    df <- df[get(recordIdColumn) %in% keepSamples, ]


    ## Format data for the different differential abundance methods.

    # First, remove id columns and any columns that are all 0s.
    cleanedData <- purrr::discard(df[, -..allIdColumns], function(col) {identical(union(unique(col), c(0, NA)), c(0, NA))})
    droppedColumns <- setdiff(names(df[, -..allIdColumns]), names(cleanedData))

    # Next, transpose abundance data to get a counts matrix with taxa as rows and samples as columns
    counts <- data.table::transpose(cleanedData)
    rownames(counts) <- names(cleanedData)
    colnames(counts) <- df[[recordIdColumn]]

    # Then, format metadata. Recall samples are rows and variables are columns
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
                                              design = as.formula(paste0("~",comparatorColName)),
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
                  fix_formula = comparatorColName, rand_formula = NULL,
                  p_adj_method = "holm", prv_cut=0,
                  group = comparatorColName)

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
    result@parameters <- paste0('recordIdColumn = ', recordIdColumn,", comparatorColName = ", comparatorColName, ', method = ', method, ', groupA =', getGroupLabels(comparator, "groupA"), ', groupB = ', getGroupLabels(comparator, "groupB"))
    result@droppedColumns <- droppedColumns


    # The resulting data should contain only the samples actually used.
    result@data <- df[, ..allIdColumns]
    names(result@data) <- stripEntityIdFromColumnHeader(names(result@data))


    validObject(result)
    veupathUtils::logWithTime(paste('Differential abundance computation completed with parameters recordIdColumn = ', recordIdColumn,", comparatorColName = ", comparatorColName, ', method = ', method, ', groupA =', getGroupLabels(comparator, "groupA"), ', groupB = ', getGroupLabels(comparator, "groupB")), verbose)
    
    return(result)
})