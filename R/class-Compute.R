# TODO tests

# so S4 will recognize data.table class as inheriting from data.frame
setOldClass(c("data.table", "data.frame"))

check_compute_result <- function(object) {
    errors <- character()
    variables <- object@computedVariableMetadata

    if (length(object@name) != 1) {
      msg <- "Compute result name must have a single value."
      errors <- c(errors, msg) 
    }
   
    # TODO think this isnt working for collections
    #col_names <- unlist(lapply(as.list(variables), function(x) {veupathUtils::getColName(x@variableSpec)}))
    #if (!all(col_names %in% names(object@data))) {
    #  msg <- paste("Some specified computed variables are not present in compute result data.frame")
    #  errors <- c(errors, msg)
    #}

    # TODO think this isnt working for the empty placeholders
    #var_classes <- unlist(lapply(as.list(variables), function(x) {x@variableClass}))
    #if (!all(var_classes %in% 'computed')) {
    #  msg <- paste("Some specified computed variables have the wrong variable class.")
    #  errors <- c(errors, msg) 
    #}

    #TODO check that data col headers dont have entityID prepended
    #TODO check col order for ids matches the order of ancestorIds, recordId

    # think we always want a data.table by now, not sure how to enforce that in the class def
    if (!'data.table' %in% class(object@data)) {
      msg <- paste("Compute result data object should be a data.table")
      errors <- c(errors, msg)
    }

    return(if (length(errors) == 0) TRUE else errors)
}

#' Compute Result
#' 
#' A class for consistently representing the results of computes. 
#' This includes their representation in R, as JSON and how they are written to files.
#' 
#' @slot data A data.frame of values where computed variables are columns and samples rows.
#' @slot recordIdColumn The name of the column containing IDs for the samples. All other columns will be treated as computed values.
#' @slot computedVariableMetadata veupathUtils::VariableMetadataList detailing the computed variables.
#' @slot computationDetails An optional message about the computed results.
#' @slot parameters A record of the input parameters used to generate the computed results.
#' 
#' @name ComputeResult-class
#' @rdname ComputeResult-class
#' @export
ComputeResult <- setClass("ComputeResult", representation(
    name = 'character',
    data = 'data.frame',
    recordIdColumn = 'character',
    ancestorIdColumns = 'character',
    computedVariableMetadata = 'VariableMetadataList',
    computationDetails = 'character',
    parameters = 'character'
), prototype = prototype(
    name = NA_character_,
    recordIdColumn = NA_character_,
    computationDetails = NA_character_,
    parameters = NA_character_
), validity = check_compute_result)

check_abundance_data <- function(object) {
    errors <- character()
    df <- object@data
    record_id_col <- object@recordIdColumn
    ancestor_id_cols <- object@ancestorIdColumns
    
    if (length(record_id_col) != 1) {
      msg <- "Record ID column must have a single value."
      errors <- c(errors, msg) 
    }

    if (!record_id_col %in% names(df)) {
      msg <- paste("Record ID column is not present in abundance data.frame")
      errors <- c(errors, msg)
    }

    if (!!length(ancestor_id_cols)) {
      if (!all(ancestor_id_cols %in% names(df))) {
        msg <- paste("Not all ancestor ID columns are present in abundance data.frame")
        errors <- c(errors, msg)
      }
    }

    #if (!all(unlist(lapply(df[, !(names(df) %in% record_id_col)], is.numeric)))) {
    #  msg <- paste("All columns except the record ID column must be numeric")
    #  errors <- c(errors, msg)
    #}

    #TODO exclude ancestors from this
    #if (uniqueN(veupathUtils::strSplit(names(df), ".", ncol=2, index=1)) > 1) {
    #  msg <- paste("All columns must belong to the same entity.")
    #  errors <- c(errors, msg)
    #}

    return(if (length(errors) == 0) TRUE else errors)
}

#' Abundance Data
#' 
#' A class for working with microbiome or ecological abundance data.
#' 
#' @slot data A data.frame of abundance values with species as columns and samples as rows
#' @slot recordIdColumn The name of the column containing IDs for the samples. All other columns will be treated as abundance values.
#' @slot imputeZero A logical indicating whether NA/ null values should be replaced with zeros.
#' 
#' @name AbundanceData-class
#' @rdname AbundanceData-class
#' @export 
AbundanceData <- setClass("AbundanceData", representation(
    data = 'data.frame',
    recordIdColumn = 'character',
    ancestorIdColumns = 'character',
    imputeZero = 'logical'
), prototype = prototype(
    recordIdColumn = NA_character_,
    imputeZero = TRUE
), validity = check_abundance_data)