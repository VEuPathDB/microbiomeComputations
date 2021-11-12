# Tests for alpha diversity functions
test_that('alphaDiv returns a correctly formatted data.table', {

  df <- testOTU
  
  results <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(typeof(results$entity.SampleID), typeof(results$entity.alphaDiversity)), c('character','double'))
  
  results <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(typeof(results$entity.SampleID), typeof(results$entity.alphaDiversity)), c('character','double'))
  
  results <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(typeof(results$entity.SampleID), typeof(results$entity.alphaDiversity)), c('character','double'))
  
})

test_that("alphaDiv returns a data.table with the correct attributes", {
  
  df <- testOTU
  
  results <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('variableId','entityId','dataType','dataShape','displayLabel','defaultRange'))
  expect_equal(attr$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Shannon Diversity')
  
  results <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('variableId','entityId','dataType','dataShape','displayLabel','defaultRange'))
  expect_equal(attr$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Simpson Diversity')
  
  results <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('variableId','entityId','dataType','dataShape','displayLabel','defaultRange'))
  expect_equal(attr$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Pielou\'s Evenness')
  
})

test_that("alphaDivApp produces an appropriately structured list of computations", {
  
  df <- testOTU
  
  # Default - use all methods
  appResults <- alphaDivApp(df, "entity.SampleID", verbose=F)
  expect_equal(length(appResults), 3)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, nrow))), nrow(df))
  
  # Test using a subset of methods
  appResults <- alphaDivApp(df, "entity.SampleID", methods=c('evenness','simpson'), verbose=F)
  expect_equal(length(appResults), 2)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, nrow))), nrow(df))
  
})

test_that("alphaDivApp output is correctly represented in json", {
  
  df <- testOTU
  
  # Default - use all methods
  nMethods <- 3
  appResults <- alphaDivApp(df, "entity.SampleID", verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('variableId','entityId','dataType','dataShape','displayLabel','defaultRange','values'))
  expect_equal(nrow(jsonList$computations$computedVariableDetails), nMethods)
  expect_equal(jsonList$computations$computedVariableDetails$variableId[[1]], 'alphaDiversity')
  expect_equal(jsonList$computations$computedVariableDetails$entityId[[1]], 'entity')
  expect_equal(jsonList$computations$computedVariableDetails$displayLabel[[1]], 'Shannon Diversity')
  expect_equal(jsonList$computations$computedVariableDetails$defaultRange[[1]], c(0, 1))
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  expect_equal(nrow(jsonList$computations$recordVariableDetails), nMethods)
  expect_equal(names(jsonList$computations$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$computations$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$computations$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$computations$recordVariableDetails$values[[1]]), nrow(df))
  
  

  
  # Test using a subset of methods
  nMethods <- 2
  appResults <- alphaDivApp(df, "entity.SampleID", methods=c('evenness','simpson'), verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('variableId','entityId','dataType','dataShape','displayLabel','defaultRange','values'))
  expect_equal(nrow(jsonList$computations$computedVariableDetails), nMethods)
  expect_equal(jsonList$computations$computedVariableDetails$variableId[[1]], 'alphaDiversity')
  expect_equal(jsonList$computations$computedVariableDetails$entityId[[1]], 'entity')
  expect_equal(jsonList$computations$computedVariableDetails$displayLabel[[1]], 'Pielou\'s Evenness')
  expect_equal(jsonList$computations$computedVariableDetails$defaultRange[[1]], c(0, 1))
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  expect_equal(nrow(jsonList$computations$recordVariableDetails), nMethods)
  expect_equal(names(jsonList$computations$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$computations$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$computations$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$computations$recordVariableDetails$values[[1]]), nrow(df))
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
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), 'entity.SampleID')
  expect_equal(typeof(results$entity.SampleID), c('character'))
  attr <- attributes(results)
  expect_equal(attr$computationDetails, "Error: alpha diversity evenness failed: input data must be non-negative")
  expect_equal(typeof(attr$parameters), 'character')
  expect_equal(typeof(attr$recordVariable), 'character')
  expect_equal(typeof(attr$computedVariableDetails$variableId), 'character')
  expect_equal(typeof(attr$computedVariableDetails$entityId), 'character')
  expect_equal(typeof(attr$computedVariableDetails$dataType), 'character')
  expect_equal(typeof(attr$computedVariableDetails$dataShape), 'character')
  expect_equal(typeof(attr$computedVariableDetails$displayLabel), 'character')
  expect_equal(typeof(attr$computedVariableDetails$defaultRange), 'double')
  
  
})