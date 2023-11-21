# considered making this an as.data.table method, but i like being consistent with the other classes+methods in this package
# plus, that wouldnt be consistent w data.table::as.data.table either really
#'@export 
setMethod("getSampleMetadata", signature("SampleMetadata"), function(object, asCopy = c(TRUE, FALSE), includeIds = c(TRUE, FALSE)) {
  asCopy <- veupathUtils::matchArg(asCopy)
  includeIds <- veupathUtils::matchArg(includeIds)
  dt <- object@data

  # Check that incoming dt meets requirements
  if (!inherits(dt, 'data.table')) {
    data.table::setDT(dt)
  }

  if (asCopy) {
    dt <- data.table::copy(dt)
  }

  if (!includeIds) {
    allIdColumns <- c(object@recordIdColumn, object@ancestorIdColumns)
    dt <- dt[, -..allIdColumns]
  }

  return(dt)
})