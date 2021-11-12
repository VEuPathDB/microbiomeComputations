# Tests for alpha diversity functions
test_that('alphaDiv returns a correctly formatted data.table', {

  df <- testOTU
  
  results <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(class(results$entity.SampleID), class(results$entity.alphaDiversity)), c('character','numeric'))
  
  results <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(class(results$entity.SampleID), class(results$entity.alphaDiversity)), c('character','numeric'))
  
  results <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(class(results$entity.SampleID), class(results$entity.alphaDiversity)), c('character','numeric'))
  
})

test_that("alphaDiv returns a data.table with the correct attributes", {
  
  df <- testOTU
  
  results <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariables') %in% names(attr)))
  expect_equal(attr$parameters, 'shannon')
  expect_equal(names(attr$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableMetadata), c('displayLabel','defaultRange'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel, 'Shannon Diversity')
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$defaultRange, c(0, 1))
  
  results <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariables') %in% names(attr)))
  expect_equal(attr$parameters, 'simpson')
  expect_equal(names(attr$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableMetadata), c('displayLabel','defaultRange'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel, 'Simpson Diversity')
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$defaultRange, c(0, 1))
  
  results <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariables') %in% names(attr)))
  expect_equal(attr$parameters, 'evenness')
  expect_equal(names(attr$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableMetadata), c('displayLabel','defaultRange'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel, 'Pielou\'s Evenness')
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$defaultRange, c(0, 1))
  
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
  appResults <- alphaDivApp(df, 'entity.SampleID', verbose=F)
  outJson <- getAppJson(appResults, 'entity.SampleID')
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets','recordVariableDetails'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  # recordVariableDetails
  expect_equal(names(jsonList$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$recordVariableDetails$values), nrow(df))
  # computations
  expect_equal(names(jsonList$computations), c('computationDetails', 'computedVariables'))
  expect_equal(nrow(jsonList$computations), nMethods)
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariables
  expect_equal(length(jsonList$computations$computedVariables), nMethods)
  expect_equal(names(jsonList$computations$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape','values'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$variableId[[1]], 'alphaDiversity')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$entityId[[1]], 'entity')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataType[[1]], 'NUMBER')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataShape[[1]], 'CONTINUOUS')
  expect_equal(length(jsonList$computations$computedVariables[[1]]$computedVariableDetails$values[[1]]), nrow(df))
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableMetadata), c('displayLabel','defaultRange'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableMetadata$displayLabel[[1]], 'Shannon Diversity')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableMetadata$defaultRange[[1]], c(0, 1))
  
  

  
  # Test using a subset of methods
  nMethods <- 2
  appResults <- alphaDivApp(df, 'entity.SampleID', methods=c('evenness','simpson'), verbose=F)
  outJson <- getAppJson(appResults, 'entity.SampleID')
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets','recordVariableDetails'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  # recordVariableDetails
  expect_equal(names(jsonList$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$recordVariableDetails$values), nrow(df))
  # computations
  expect_equal(names(jsonList$computations), c('computationDetails', 'computedVariables'))
  expect_equal(nrow(jsonList$computations), nMethods)
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariables
  expect_equal(length(jsonList$computations$computedVariables), nMethods)
  expect_equal(names(jsonList$computations$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape','values'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$variableId[[1]], 'alphaDiversity')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$entityId[[1]], 'entity')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataType[[1]], 'NUMBER')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataShape[[1]], 'CONTINUOUS')
  expect_equal(length(jsonList$computations$computedVariables[[1]]$computedVariableDetails$values[[1]]), nrow(df))
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableMetadata), c('displayLabel','defaultRange'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableMetadata$displayLabel[[1]], 'Pielou\'s Evenness')
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableMetadata$defaultRange[[1]], c(0, 1))
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
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$variableId), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$entityId), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$dataType), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$dataShape), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableMetadata$defaultRange), 'double')
  
  
})