# Tests for beta diversity functions
test_that('betaDiv returns a correctly formatted data.table', {
  
  df <- testOTU

  results <- betaDiv(df, "entity.SampleID", method='bray', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(results, class))), c('character','numeric','numeric'))
  
  results <- betaDiv(df, "entity.SampleID", method='jaccard', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(results, class))), c('character','numeric','numeric'))
  
  results <- betaDiv(df, "entity.SampleID", method='jsd', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(results, class))), c('character','numeric','numeric'))
  
  
})

test_that("betaDiv returns a data.table with the correct attributes" , {
  
  df <- testOTU
  
  results <- betaDiv(df, "entity.SampleID", method='bray', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariables','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'bray')
  expect_equal(names(attr$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableMetadata), c('displayLabel'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$variableId, c('Axis1','Axis2'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$entityId, rep('entity',2))
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel, c('Axis1 15.3%','Axis2 5.7%'))
  
  results <- betaDiv(df, "entity.SampleID", method='jaccard', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariables','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'jaccard')
  expect_equal(names(attr$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableMetadata), c('displayLabel'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$variableId, c('Axis1','Axis2'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$entityId, rep('entity',2))
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel, c('Axis1 10.0%','Axis2 4.3%'))
  
  results <- betaDiv(df, "entity.SampleID", method='jsd', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariables','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'jsd')
  expect_equal(names(attr$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariables[[1]]$computedVariableMetadata), c('displayLabel'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$variableId, c('Axis1','Axis2'))
  expect_equal(attr$computedVariables[[1]]$computedVariableDetails$entityId, rep('entity',2))
  expect_equal(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel,c('Axis1 25.2%','Axis2 17.5%'))

})


test_that("betaDivApp produces an appropriately structured list of computations", {
  
  df <- testOTU
  
  # Default - use all methods
  appResults <- betaDivApp(df, "entity.SampleID", verbose=F)
  expect_equal(length(appResults), 3)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, nrow))), nrow(df))
  
  # Test using a subset of methods
  appResults <- betaDivApp(df, "entity.SampleID", methods = c('jaccard','bray'), verbose=F)
  expect_equal(length(appResults), 2)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, nrow))), nrow(df))
  
})


test_that("betaDivApp output is correctly represented in json", {
  
  df <- testOTU
  
  nMethods <- 1
  appResults <- betaDivApp(df, "entity.SampleID", methods=c('bray'), verbose=F)
  outJson <- getAppJson(appResults, "entity.SampleID")
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets','recordVariableDetails'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  # recordVariableDetails
  expect_equal(names(jsonList$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$recordVariableDetails$values), nrow(df))
  # computations
  expect_equal(names(jsonList$computations), c('computationDetails', 'pcoaVariance', 'computedVariables'))
  expect_equal(nrow(jsonList$computations), nMethods)
  expect_equal(jsonList$computations$pcoaVariance[[1]], c(15.3, 5.7))
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariables
  expect_equal(length(jsonList$computations$computedVariables), nMethods)
  expect_equal(names(jsonList$computations$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection','values'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$variableId[[1]], c('Axis1', 'Axis2'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$entityId[[1]], c('entity','entity'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataType[[1]], c('NUMBER', 'NUMBER'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataShape[[1]], c('CONTINUOUS', 'CONTINUOUS'))
  expect_equal(ncol(jsonList$computations$computedVariables[[1]]$computedVariableDetails$values[[1]]), nrow(df))
  expect_false(jsonList$computations$computedVariables[[1]]$computedVariableDetails$isCollection[[1]])
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableMetadata), c('displayLabel'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableMetadata$displayLabel[[1]], c('Axis1 15.3%','Axis2 5.7%'))
  
  
  nMethods <- 2
  appResults <- betaDivApp(df, "entity.SampleID", methods = c('bray','jaccard'), verbose=F)
  outJson <- getAppJson(appResults, "entity.SampleID")
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets','recordVariableDetails'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  # recordVariableDetails
  expect_equal(names(jsonList$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$recordVariableDetails$values), nrow(df))
  # computations
  expect_equal(names(jsonList$computations), c('computationDetails', 'pcoaVariance', 'computedVariables'))
  expect_equal(nrow(jsonList$computations), nMethods)
  expect_equal(jsonList$computations$pcoaVariance[[1]], c(15.3, 5.7))
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariables
  expect_equal(length(jsonList$computations$computedVariables), nMethods)
  expect_equal(names(jsonList$computations$computedVariables[[1]]), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection','values'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$variableId[[1]], c('Axis1', 'Axis2'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$entityId[[1]], c('entity','entity'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataType[[1]], c('NUMBER', 'NUMBER'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableDetails$dataShape[[1]], c('CONTINUOUS', 'CONTINUOUS'))
  expect_equal(ncol(jsonList$computations$computedVariables[[1]]$computedVariableDetails$values[[1]]), nrow(df))
  expect_false(jsonList$computations$computedVariables[[1]]$computedVariableDetails$isCollection[[1]])
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariables[[1]]$computedVariableMetadata), c('displayLabel'))
  expect_equal(jsonList$computations$computedVariables[[1]]$computedVariableMetadata$displayLabel[[1]], c('Axis1 15.3%','Axis2 5.7%'))
  
  
})

test_that("betaDiv results are consistent", {
  
  expect_snapshot_value({
    df <- testOTU
    appResults <- betaDivApp(df, "entity.SampleID", verbose=F)
  }, style = "serialize")
  
})

test_that("betaDiv fails gracefully", {
  
  df <- testOTU
  df$entity.Abiotrophia <- NA
  
  results <- betaDiv(df, "entity.SampleID", method='bray', verbose=T)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), 'entity.SampleID')
  expect_equal(typeof(results$entity.SampleID), c('character'))
  attr <- attributes(results)
  expect_equal(attr$computationDetails, "Error: beta diversity bray failed: missing values are not allowed with argument 'na.rm = FALSE'")
  expect_equal(typeof(attr$parameters), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$variableId), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$entityId), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$dataType), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableDetails$dataShape), 'character')
  expect_equal(typeof(attr$computedVariables[[1]]$computedVariableMetadata$displayLabel), 'character')
  expect_equal(typeof(attr$pcoaVariance), 'double')
  
})