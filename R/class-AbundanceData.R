check_abundance_data <- function(object) {
    errors <- character()
    df <- object@data
    record_id_col <- object@recordIdColumn
    ancestor_id_cols <- object@ancestorIdColumns

    # Ensure id columns are valid
    msg <- validateIdColumns(df, record_id_col, ancestor_id_cols)
    errors <- c(errors, msg)

    # Abundance data should all come from the same entity
    # using the presence of the period to indicate eda services formatted data
    if (all(grepl(".", names(df), fixed = TRUE))) {
      if (uniqueN(veupathUtils::strSplit(names(df)[!names(df) %in% ancestor_id_cols], ".", ncol=2, index=1)) > 1) {
        msg <- paste("All columns must belong to the same entity.")
        errors <- c(errors, msg)
      }
    }

    allDataColsNumeric <- all(unlist(lapply(df[, !(names(df) %in% c(record_id_col, ancestor_id_cols))], is.numeric)))
    if (inherits(df, 'data.table')) allDataColsNumeric <- all(unlist(lapply(df[, !(names(df) %in% c(record_id_col, ancestor_id_cols)), with=F], is.numeric)))
    if (!allDataColsNumeric) {
      msg <- paste("All columns except the ID columns must be numeric.")
      errors <- c(errors, msg)
    }

    if (any(df < 0, na.rm=TRUE)) {
      msg <- paste("Abundance data cannot contain negative values.")
      errors <- c(errors, msg)
    }

    if (!!length(object@sampleMetadata@data)) {
      sampleMetadata <- object@sampleMetadata

      if (!identical(sampleMetadata@recordIdColumn, record_id_col)) {
        msg <- paste("Records in the sample metadata and the abundance data refer to different entities.")
      }
      if (!setequal(sampleMetadata@data[[sampleMetadata@recordIdColumn]], df[[record_id_col]])) {
        msg <- paste("Samples do not match between the sample metadata and abundance data.")
        errors <- c(errors, msg)
      }
      if (!identical(sampleMetadata@data[[sampleMetadata@recordIdColumn]], df[[record_id_col]])) {
        msg <- paste("Samples in the sample metadata are not in the same order as in the abundance data.")
        errors <- c(errors, msg)
      }
      if (setequal(names(sampleMetadata@data), c(record_id_col, ancestor_id_cols))) {
        msg <- paste("The sample metadata only contains record ID and ancestor ID columns but no metadata variables.")
        errors <- c(errors, msg)
      }
    }

    return(if (length(errors) == 0) TRUE else errors)
}

#' Abundance Data
#' 
#' A class for working with microbiome or ecological abundance data.
#' 
#' @slot data A data.frame of abundance values with species as columns and samples as rows
#' @slot sampleMetadata A SampleMetadata object of metadata about the samples with samples as rows and metadata variables as columns
#' @slot recordIdColumn The name of the column containing IDs for the samples. All other columns will be treated as abundance values.
#' @slot ancestorIdColumns A character vector of column names representing parent entities of the recordIdColumn.
#' @slot imputeZero A logical indicating whether NA/ null values should be replaced with zeros.
#' @slot removeEmptySamples A logical indicating whether empty (all NA/ zero)samples should be removed.
#' @name AbundanceData-class
#' @rdname AbundanceData-class
#' @include class-SampleMetadata.R
#' @export 
AbundanceData <- setClass("AbundanceData", representation(
    data = 'data.frame',
    sampleMetadata = 'SampleMetadata',
    recordIdColumn = 'character',
    ancestorIdColumns = 'character',
    imputeZero = 'logical',
    removeEmptySamples = 'logical'
), prototype = prototype(
    recordIdColumn = NA_character_,
    imputeZero = TRUE,
    removeEmptySamples = TRUE
), validity = check_abundance_data)