check_sample_metadata <- function(object) {
    errors <- character()
    df <- object@data
    record_id_col <- object@recordIdColumn
    ancestor_id_cols <- object@ancestorIdColumns
    
    if (length(record_id_col) != 1) {
      msg <- "Record ID column must have a single value."
      errors <- c(errors, msg) 
    }

    if (!record_id_col %in% names(df)) {
      msg <- paste("Record ID column is not present in sample metadata data.frame")
      errors <- c(errors, msg)
    }

    if (!!length(ancestor_id_cols)) {
      if (!all(ancestor_id_cols %in% names(df))) {
        msg <- paste("Not all ancestor ID columns are present in sample metadata data.frame")
        errors <- c(errors, msg)
      }
    }
    

    return(if (length(errors) == 0) TRUE else errors)
}

#' Sample Metadata
#' 
#' A class for working with metadata.
#' 
#' @slot data A data.frame of metadata variables with samples as rows and variables as columns
#' @slot sampleMetadata A data.frame of metadata about the samples with samples as rows and metadata variables as columns
#' @slot recordIdColumn The name of the column containing IDs for the samples. 
#' @slot ancestorIdColumns A character vector of column names representing parent entities of the recordIdColumn.
#' @slot variableMetadata A VariableMetadataList containing information about the variables contained in data.
#' @name SampleMetadata-class
#' @rdname SampleMetadata-class
#' @export 
SampleMetadata <- setClass("SampleMetadata", representation(
    data = 'data.frame',
    recordIdColumn = 'character',
    ancestorIdColumns = 'character',
    variableMetadata = 'VariableMetadataList'
), prototype = prototype(
    recordIdColumn = NA_character_
), validity = check_sample_metadata)