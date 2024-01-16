#' Get data.table of abundances from AbundanceData
#'
#' Returns a data.table of abundances, respecting the
#' `imputeZero` slot.
#' 
#' @param object AbundanceData
#' @param ignoreImputeZero boolean indicating whether we should respect the imputeZero slot
#' @param includeIds boolean indicating whether we should include recordIdColumn and ancestorIdColumns
#' @return data.table of abundances
#' @export
setGeneric("getAbundances",
  function(object, ignoreImputeZero = c(FALSE, TRUE), includeIds = c(TRUE, FALSE)) standardGeneric("getAbundances"),
  signature = c("object")
)

#'@export 
setMethod("getAbundances", signature("AbundanceData"), function(object, ignoreImputeZero = c(FALSE, TRUE), includeIds = c(TRUE, FALSE)) {
  ignoreImputeZero <- veupathUtils::matchArg(ignoreImputeZero)
  includeIds <- veupathUtils::matchArg(includeIds)
  dt <- object@data
  allIdColumns <- c(object@recordIdColumn, object@ancestorIdColumns)

  # Check that incoming dt meets requirements
  if (!inherits(dt, 'data.table')) {
    # this might technically be bad form, but i think its ok in this context
    data.table::setDT(dt)
  }

  if (object@removeEmptySamples) {
    dt.noIds <- dt[, -..allIdColumns]
    # Remove samples with NA in all columns
    dt <- dt[rowSums(is.na(dt.noIds)) != ncol(dt.noIds),]

    # Remove samples with zero in all columns
    dt <- dt[rowSums(dt.noIds == 0) != ncol(dt.noIds),]
  }

  # Replace NA values with 0
  if (!ignoreImputeZero && object@imputeZero) {
    veupathUtils::setNaToZero(dt)
  }

  if (!includeIds) {
    dt <- dt[, -..allIdColumns]
  }

  return(dt)
})

#' Get data.table of sample metadata from AbundanceData
#'
#' Returns a data.table of sample metadata
#' 
#' @param object AbundanceData
#' @param asCopy boolean indicating whether to return the data as a copy or by reference
#' @param includeIds boolean indicating whether we should include recordIdColumn and ancestorIdColumns
#' @return data.table of sample metadata
#' @export
setGeneric("getSampleMetadata",
  function(object, asCopy = c(TRUE, FALSE), includeIds = c(TRUE, FALSE)) standardGeneric("getSampleMetadata"),
  signature = c("object")
)

#'@export 
setMethod("getSampleMetadata", signature("AbundanceData"), function(object, asCopy = c(TRUE, FALSE), includeIds = c(TRUE, FALSE)) {
  asCopy <- veupathUtils::matchArg(asCopy)
  includeIds <- veupathUtils::matchArg(includeIds)
  dt <- object@sampleMetadata@data
  allIdColumns <- c(object@recordIdColumn, object@ancestorIdColumns)

  # Check that incoming dt meets requirements
  if (!inherits(dt, 'data.table')) {
    data.table::setDT(dt)
  }

  if (asCopy) {
    dt <- data.table::copy(dt)
  }

  if (object@removeEmptySamples) {
    # not using getAbundances here bc i want the empty samples here
    abundances <- object@data[, -..allIdColumns]

    # Remove metadata for samples with NA in all columns
    dt <- dt[rowSums(is.na(abundances)) != ncol(abundances),]

    # Remove metadata for samples with zero in all columns
    dt <- dt[rowSums(abundances == 0) != ncol(abundances),]
  }

  if (!includeIds) {
    dt <- dt[, -..allIdColumns]
  }

  return(dt)
})

#' Drop samples with incomplete SampleMetadata
#'
#' Modifies the data and sampleMetadata slots of an 
#' AbundanceData object, to exclude samples with 
#' missing SampleMetadata for a specified column.
#' 
#' @param object AbundanceData
#' @param colName String providing the column name in SampleMetadata to check for completeness
#' @param verbose boolean indicating if timed logging is desired
#' @return AbundanceData with modified data and sampleMetadata slots
#' @export
setGeneric("removeIncompleteSamples",
  function(object, colName = character(), verbose = c(TRUE, FALSE)) standardGeneric("removeIncompleteSamples"),
  signature = c("object")
)

#'@export 
setMethod("removeIncompleteSamples", signature("AbundanceData"), function(object, colName = character(), verbose = c(TRUE, FALSE)) {
  df <- getAbundances(object)
  sampleMetadata <- getSampleMetadata(object)

  # Remove samples with NA from data and metadata
  if (any(is.na(sampleMetadata[[colName]]))) {
    veupathUtils::logWithTime("Found NAs in specified variable. Removing these samples.", verbose)
    samplesWithData <- which(!is.na(sampleMetadata[[colName]]))
    # Keep samples with data. Recall the AbundanceData object requires samples to be in the same order
    # in both the data and metadata
    sampleMetadata <- sampleMetadata[samplesWithData, ]
    df <- df[samplesWithData, ]

    object@data <- df
    object@sampleMetadata <- SampleMetadata(
      data = sampleMetadata,
      recordIdColumn = object@sampleMetadata@recordIdColumn
    )
    validObject(object)
  }

  return(object)
})