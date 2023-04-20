check_absolute_abundance_data <- function(object) {
    errors <- character()
    df <- object@data
    record_id_col <- object@recordIdColumn

    # otu <- df[, -..record_id_col]

    # if (!identical(otu, round(otu))) {
    #   msg <- "Absolute abundance data must be integer numbers."
    #   errors <- c(errors, msg)
    # }
    

    return(if (length(errors) == 0) TRUE else errors)
}

#' Absolute Abundance Data
#' 
#' A class for working with count data from a microbial or ecological assay.
#' 
#' @slot data A data.frame of integer abundance counts with species as columns and samples as rows
#' @slot sampleMetadata A data.frame of metadata about the samples with samples as rows and metadata variables as columns
#' @slot recordIdColumn The name of the column containing IDs for the samples. All other columns will be treated as abundance values.
#' @slot ancestorIdColumns A character vector of column names representing parent entities of the recordIdColumn.
#' @slot imputeZero A logical indicating whether NA/ null values should be replaced with zeros.
#' @name AbsoluteAbundanceData-class
#' @rdname AbsoluteAbundanceData-class
#' @include class-AbundanceData.R
#' @export 
AbsoluteAbundanceData <- setClass("AbsoluteAbundanceData", contains = "AbundanceData", validity = check_absolute_abundance_data)