check_sample_metadata <- function(object) {
    errors <- character()
    df <- object@data
    record_id_col <- object@recordIdColumn
    ancestor_id_cols <- object@ancestorIdColumns
    
    msg <- validateIdColumns(df, record_id_col, ancestor_id_cols)
    

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
#' @name SampleMetadata-class
#' @rdname SampleMetadata-class
#' @export 
SampleMetadata <- setClass("SampleMetadata", representation(
    data = 'data.frame',
    recordIdColumn = 'character',
    ancestorIdColumns = 'character'
), prototype = prototype(
    recordIdColumn = NA_character_
), validity = check_sample_metadata)