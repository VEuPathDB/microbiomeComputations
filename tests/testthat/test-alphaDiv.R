# Tests for alpha diversity functions
test_that('alphaDiv returns a correctly formatted data.table', {

  df <- testOTU
  
  results <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(typeof(results$entity.SampleID), typeof(results$entity.alphaDiversity)), c('character','double'))
  
  results <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(typeof(results$entity.SampleID), typeof(results$entity.alphaDiversity)), c('character','double'))
  
  results <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(typeof(results$entity.SampleID), typeof(results$entity.alphaDiversity)), c('character','double'))
  
})

test_that("alphaDiv returns a data.table with the correct attributes", {
  
  df <- testOTU
  
  results <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('id','entity','displayLabel','defaultRange'))
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$entity, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Shannon')
  
  results <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('id','entity','displayLabel','defaultRange'))
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$entity, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Simpson')
  
  results <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('id','entity','displayLabel','defaultRange'))
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$entity, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Pielou\'s Evenness')
  
})

test_that("alphaDivApp produces an appropriately structured list of computations", {
  
  df <- testOTU
  
  # Default - use all methods
  appResults <- alphaDivApp(df, "entity.SampleID", verbose=F)
  expect_equal(length(appResults), 3)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, NROW))), NROW(df))
  
  # Test using a subset of methods
  appResults <- alphaDivApp(df, "entity.SampleID", methods=c('evenness','simpson'), verbose=F)
  expect_equal(length(appResults), 2)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, NROW))), NROW(df))
  
})

test_that("alphaDivApp output is correctly represented in json", {
  
  df <- testOTU
  
  # Default - use all methods
  appResults <- alphaDivApp(df, "entity.SampleID", verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 3)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','defaultRange','values'))
  expect_equal(jsonList$computations$computedVariableDetails$id[[1]], 'alphaDiversity')

  
  # Test using a subset of methods
  appResults <- alphaDivApp(df, "entity.SampleID", methods=c('evenness','simpson'), verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 2)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','defaultRange','values'))
  expect_equal(jsonList$computations$computedVariableDetails$id[[1]], 'alphaDiversity')
})

test_that("alphaDiv results are consistent", {
  
  expect_snapshot_value({
    df <- testOTU
    appResults <- alphaDivApp(df, "entity.SampleID", verbose=F)
  }, style = "serialize")
  
})

test_that("alphaDiv fails gracefully", {
  
  df <- testOTU
  df$entity.Abiotrophia <- -1
  
  results <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=T)
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), 'entity.SampleID')
  expect_equal(typeof(results$entity.SampleID), c('character'))
  attr <- attributes(results)
  expect_equal(attr$computationDetails, "Error: alpha diversity evenness failed: input data must be non-negative")
  expect_equal(typeof(attr$parameters), 'character')
  expect_equal(typeof(attr$recordVariable), 'character')
  expect_equal(typeof(attr$computedVariableDetails$id), 'character')
  expect_equal(typeof(attr$computedVariableDetails$entity), 'character')
  expect_equal(typeof(attr$computedVariableDetails$displayLabel), 'character')
  expect_equal(typeof(attr$computedVariableDetails$defaultRange), 'double')
  
  
})