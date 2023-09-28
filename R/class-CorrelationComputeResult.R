check_correlation_compute_result <- function(object) {
    errors <- character()

# ??????
    # if (is.na(object@name)) {
    #   msg <- "Compute result must have a name."
    #   errors <- c(errors, msg)
    # } else if (length(object@name) != 1) {
    #   msg <- "Compute result name must have a single value."
    #   errors <- c(errors, msg) 
    # }

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
#' @slot data1Metadata A VariableMetadata object describing the data1 object sent to the correlation function
#' @slot data1Variables A VariableMetadataList object that describes each variable in data1.
#' @slot data2Metadata A VariableMetadata object describing the data2 object sent to the correlation function
#' @slot data2Variables A VariableMetadataList object that describes each variable in data2.
#' @name CorrelationComputeResult-class
#' @rdname CorrelationComputeResult-class
#' @export
CorrelationComputeResult <- setClass("CorrelationComputeResult", contains="ComputeResult", representation(
    data1Metadata = "VariableMetadata",
    data1Variables = "VariableMetadataList",
    data2Metadata = "VariableMetadata",
    data2Variables = "VariableMetadataList"
), prototype = prototype(
    name = NA_character_,
    recordIdColumn = NA_character_,
    computationDetails = NA_character_,
    parameters = NA_character_
), validity = check_compute_result)