# TODO roxygen docs

# so S4 will recognize data.table class as inheriting from data.frame
setOldClass(c("data.table", "data.frame"))

check_compute_result <- function(object) {
    errors <- character()
    variables <- object@computedVariableMetadata
   
    col_names <- unlist(lapply(as.list(variables), function(x) {veupathUtils::getColName(x@variableSpec)}))
    if (!all(col_names %in% names(object@data))) {
      msg <- paste("Some specified computed variables are not present in compute result data.frame")
      errors <- c(errors, msg)
    }

    var_classes <- unlist(lapply(as.list(variables), function(x) {x@variableClass}))
    if (!all(var_classes %in% 'computed')) {
      msg <- paste("Some specified computed variables have the wrong variable class.")
      errors <- c(errors, msg) 
    }

    return(if (length(errors) == 0) TRUE else errors)
}

#' @export
ComputeResult <- setClass("ComputeResult", representation(
    data = 'data.frame',
    computedVariableMetadata = 'VariableMetadataList'
), validity = check_compute_result)

check_abundance_data <- function(object) {
    errors <- character()
    record_id_col <- object@recordIdColumn
    
    if (length(record_id_col) != 1) {
      msg <- "Record ID column must have a single value."
      errors <- c(errors, msg) 
    }

    if (!record_id_col %in% names(object@data)) {
      msg <- paste("Record ID column is not present in abundance data.frame")
      errors <- c(errors, msg)
    }

    return(if (length(errors) == 0) TRUE else errors)
}

#' @export 
AbundanceData <- setClass("AbundanceData", representation(
    data = 'data.frame',
    recordIdColumn = 'character'
), prototype = prototype(
    recordIdColumn = NA_character_
), validity = check_abundance_data)