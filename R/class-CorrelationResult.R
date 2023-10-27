check_correlation_result <- function(object) {
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


#' Corrleation Result
#' 
#' A class to represent the computed results from a correlation calculation.
#' This includes their representation in R, as JSON, and how they are written to files.
#' 
#' @slot statistics An optional slot of any values. List or data.frame are recommended. It is not required to have rows or cols map to samples.
#' @slot data1Metadata A character describing data1. Should at least describe the name of data1.
#' @slot data2Metadata A character describing data2. Should at least describe the name of data2.
#' @name CorrelationResult-class
#' @rdname CorrelationResult-class
#' @export
CorrelationResult <- setClass("CorrelationResult", representation(
    data1Metadata = "character",
    data2Metadata = "character",
    statistics = "data.frame"
), prototype = prototype(
    data1Metadata = NA_character_,
    data2Metadata = NA_character_
), validity = check_correlation_result)