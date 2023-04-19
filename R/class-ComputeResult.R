# so S4 will recognize data.table class as inheriting from data.frame
setOldClass(c("data.table", "data.frame"))

check_compute_result <- function(object) {
    errors <- character()
    variables <- object@computedVariableMetadata

    if (is.na(object@name)) {
      msg <- "Compute result must have a name."
      errors <- c(errors, msg)
    } else if (length(object@name) != 1) {
      msg <- "Compute result name must have a single value."
      errors <- c(errors, msg) 
    }

    if (!length(variables)) {
      msg <- "Compute result must include computed variable metadata."
      errors <- c(errors, msg)      
    }
   
    col_names <- stripEntityIdFromColumnHeader(veupathUtils::findAllColNames(variables))
    if (!all(col_names %in% names(object@data))) {
      msg <- paste("Some specified computed variables are not present in compute result data.frame")
      errors <- c(errors, msg)
    }

    var_classes <- unlist(lapply(as.list(variables), function(x) {x@variableClass@value}))
    if (!all(var_classes %in% 'computed')) {
      msg <- paste("Some specified computed variables have the wrong variable class.")
      errors <- c(errors, msg) 
    }

    if (any(grepl(".", names(object@data), fixed = TRUE))) {
      msg <- paste("Column headers appear to be in dot notation [entityId.variableId]. They should be the raw variableId.")
      errors <- c(errors, msg)
    }
    
    expectedOutputIdColHeaders <- stripEntityIdFromColumnHeader(c(object@recordIdColumn, object@ancestorIdColumns))
    actualOutputIdColHeaders <- names(object@data)[1:length(expectedOutputIdColHeaders)]
    if (all(expectedOutputIdColHeaders != actualOutputIdColHeaders)) {
      msg <- paste("Columns must be ordered by recordIdColumn, ancestorIdColumns, and then data columns.")
      errors <- c(errors, msg) 
    }

    # think we always want a data.table by now, not sure how to enforce that in the class def
    #if (!'data.table' %in% class(object@data)) {
    #  msg <- paste("Compute result data object should be a data.table")
    #  errors <- c(errors, msg)
    #}

    return(if (length(errors) == 0) TRUE else errors)
}

#' Compute Result
#' 
#' A class for consistently representing the results of computes. 
#' This includes their representation in R, as JSON and how they are written to files.
#' 
#' @slot data A data.frame of values where computed variables are columns and samples rows.
#' @slot name The name of the compute, ex: 'alphaDiv'.
#' @slot recordIdColumn The name of the column containing IDs for the samples. All other columns will be treated as computed values.
#' @slot ancestorIdColumns A character vector of column names representing parent entities of the recordIdColumn.
#' @slot computedVariableMetadata veupathUtils::VariableMetadataList detailing the computed variables.
#' @slot statistics An optional data.frame of values. This data frame is not required to have rows or cols map to samples.
#' @slot computedStatisticMetadata An optional veupathUtils::VariableMetadataList containing details of any computed statistics
#' @slot computationDetails An optional message about the computed results.
#' @slot parameters A record of the input parameters used to generate the computed results.
#' @name ComputeResult-class
#' @rdname ComputeResult-class
#' @export
ComputeResult <- setClass("ComputeResult", representation(
    name = 'character',
    data = 'data.frame',
    recordIdColumn = 'character',
    ancestorIdColumns = 'character',
    computedVariableMetadata = 'VariableMetadataList',
    statistics = 'data.frame',
    computedStatisticMetadata = 'VariableMetadataList',
    computationDetails = 'character',
    parameters = 'character'
), prototype = prototype(
    name = NA_character_,
    recordIdColumn = NA_character_,
    computationDetails = NA_character_,
    parameters = NA_character_
), validity = check_compute_result)