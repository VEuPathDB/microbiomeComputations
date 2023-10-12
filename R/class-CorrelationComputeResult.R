check_correlation_compute_result <- function(object) {
    errors <- character()

    if (!length(object@data1Metadata) > 0) {
      msg <- "Compute result must include information about data1 in the data1Metadata slot"
      errors <- c(errors, msg)
    } else if (!length(object@data2Metadata) > 0) {
      msg <- "Compute result must include information about data2 in the data2Metadata slot"
      errors <- c(errors, msg) 
    }

    return(if (length(errors) == 0) TRUE else errors)
}


#' Corrleation Compute Result
#' 
#' A class to represent the results of a correlation calculation. Contains ComputeResult.
#' This includes their representation in R, as JSON and how they are written to files.
#' 
#' @slot data A data.frame of values where computed variables are columns and samples rows.
#' @slot name The name of the compute, ex: 'alphaDiv'.
#' @slot recordIdColumn The name of the column containing IDs for the samples. All other columns will be treated as computed values.
#' @slot ancestorIdColumns A character vector of column names representing parent entities of the recordIdColumn.
#' @slot computedVariableMetadata veupathUtils::VariableMetadataList detailing the computed variables.
#' @slot statistics An optional slot of any values. List or data.frame are recommended. It is not required to have rows or cols map to samples.
#' @slot computationDetails An optional message about the computed results.
#' @slot parameters A record of the input parameters used to generate the computed results.
#' @slot droppedColumns A character vector of column names that have been dropped for being unsuitable for the computation.
#' @slot data1Metadata A character describing data1. Should at least describe the name of data1.
#' @slot data2Metadata A character describing data2. Should at least describe the name of data2.
#' @name CorrelationComputeResult-class
#' @rdname CorrelationComputeResult-class
#' @export
CorrelationComputeResult <- setClass("CorrelationComputeResult", contains="ComputeResult", representation(
    data1Metadata = "character",
    data2Metadata = "character"
), prototype = prototype(
    name = NA_character_,
    recordIdColumn = NA_character_,
    computationDetails = NA_character_,
    parameters = NA_character_,
    data1Metadata = NA_character_,
    data2Metadata = NA_character_
), validity = check_compute_result)