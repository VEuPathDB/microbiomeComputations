# Tests for beta diversity functions
test_that('betaDiv returns a correctly formatted data.table', {
  
  df <- testOTU

  dt <- betaDiv(df, "entity.SampleID", method='bray', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(dt, class))), c('character','numeric','numeric'))
  
  dt <- betaDiv(df, "entity.SampleID", method='jaccard', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(dt, class))), c('character','numeric','numeric'))
  
  dt <- betaDiv(df, "entity.SampleID", method='jsd', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(dt, class))), c('character','numeric','numeric'))
  
  # With NAs
  nNAs <- 20
  df[sample(1:nrow(df), size=nNAs, replace = F), 2] <- NA
  expect_error(betaDiv(df, "entity.SampleID", method='jsd', naToZero=F, verbose=F))  # all three methods err
  dt <- betaDiv(df, "entity.SampleID", method='bray', verbose=F) 
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(dt, class))), c('character','numeric','numeric'))
  expect_true(sum(dt > 0) > 100)   # ensure output is not all 0
  
})

test_that("betaDiv returns a data.table with the correct attributes" , {
  
  df <- testOTU
  
  results <- betaDiv(df, "entity.SampleID", method='bray', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'bray')
  expect_equal(inherits(attr$computedVariable, "VariableMetadataList"), TRUE)
  expect_equal(attr$computedVariable[[1]]@variableSpec@variableId, 'Axis1')
  expect_equal(attr$computedVariable[[1]]@variableSpec@entityId, 'entity')
  expect_equal(attr$computedVariable[[1]]@displayName, 'Axis1 15.3%')
  
  results <- betaDiv(df, "entity.SampleID", method='jaccard', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'jaccard')
  expect_equal(inherits(attr$computedVariable, "VariableMetadataList"), TRUE)
  expect_equal(attr$computedVariable[[1]]@variableSpec@variableId, 'Axis1')
  expect_equal(attr$computedVariable[[1]]@variableSpec@entityId, 'entity')
  expect_equal(attr$computedVariable[[1]]@displayName, 'Axis1 10.0%')
  
  results <- betaDiv(df, "entity.SampleID", method='jsd', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'jsd')
  expect_equal(inherits(attr$computedVariable, "VariableMetadataList"), TRUE)
  expect_equal(attr$computedVariable[[1]]@variableSpec@variableId, 'Axis1')
  expect_equal(attr$computedVariable[[1]]@variableSpec@entityId, 'entity')
  expect_equal(attr$computedVariable[[1]]@displayName, 'Axis1 25.2%')

})



# now done in veupathUtils
#test_that("getMetadata() output is correctly represented in json", {
#  
#  df <- testOTU
#  
#  nMethods <- 1
#  results <- betaDiv(df, "entity.SampleID", method=c('bray'), verbose=F)
#  outJson <- getMetadata(results)
#  jsonList <- jsonlite::fromJSON(outJson)
#  expect_equal(names(jsonList), c('computedVariableDetails','computedVariableMetadata'))
#  # computedVariableDetails
#  expect_equal(names(jsonList$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
#  expect_equal(jsonList@variableSpec@variableId, c('Axis1', 'Axis2'))
#  expect_equal(jsonList@variableSpec@entityId, c('entity','entity'))
#  expect_equal(jsonList@variableSpec@dataType, c('NUMBER', 'NUMBER'))
#  expect_equal(jsonList@variableSpec@dataShape, c('CONTINUOUS', 'CONTINUOUS'))
#  expect_false(jsonList@variableSpec@isCollection)
#  # computedVariableMetadata
#  expect_equal(names(jsonList$computedVariableMetadata), c('displayName'))
#  expect_equal(jsonList@displayName, c('Axis1 15.3%','Axis2 5.7%'))
#
#})


test_that("betaDiv fails gracefully", {
  
  df <- testOTU
  df$entity.Abiotrophia <- NA
  
  dt <- betaDiv(df, "entity.SampleID", method='bray', naToZero=F, verbose=T)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), 'entity.SampleID')
  expect_equal(typeof(dt$entity.SampleID), c('character'))
  attr <- attributes(dt)
  expect_equal(attr$computationDetails, "Error: beta diversity bray failed: missing values are not allowed with argument 'na.rm = FALSE'")
  expect_equal(typeof(attr$parameters), 'character')
  expect_equal(typeof(attr$computedVariable@variableSpec@variableId), 'character')
  expect_equal(typeof(attr$computedVariable@variableSpec@entityId), 'character')
  expect_equal(typeof(attr$computedVariable@dataType@value), 'character')
  expect_equal(typeof(attr$computedVariable@dataShape@value), 'character')
  expect_equal(typeof(attr$computedVariable@displayName), 'character')
  expect_equal(typeof(attr$pcoaVariance), 'double')
  
})