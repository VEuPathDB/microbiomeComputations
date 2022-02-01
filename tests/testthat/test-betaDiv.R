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
  expect_error(betaDiv(df, "entity.SampleID", method='jsd', verbose=F))  # all three methods err
  dt <- betaDiv(df, "entity.SampleID", method='bray', naToZero=T, verbose=F) 
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
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('displayName'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, c('Axis1','Axis2'))
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, rep('entity',2))
  expect_equal(attr$computedVariable$computedVariableMetadata$displayName, c('Axis1 15.3%','Axis2 5.7%'))
  
  results <- betaDiv(df, "entity.SampleID", method='jaccard', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'jaccard')
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('displayName'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, c('Axis1','Axis2'))
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, rep('entity',2))
  expect_equal(attr$computedVariable$computedVariableMetadata$displayName, c('Axis1 10.0%','Axis2 4.3%'))
  
  results <- betaDiv(df, "entity.SampleID", method='jsd', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','pcoaVariance') %in% names(attr)))
  expect_equal(attr$parameters, 'jsd')
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('displayName'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, c('Axis1','Axis2'))
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, rep('entity',2))
  expect_equal(attr$computedVariable$computedVariableMetadata$displayName,c('Axis1 25.2%','Axis2 17.5%'))

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

  # With NAs
  nNAs <- 20
  df[sample(1:nrow(df), size=nNAs, replace = F), 2] <- NA
  appResults <- betaDivApp(df, "entity.SampleID", methods = c('jaccard','bray'), naToZero = T, verbose=F)
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
  expect_equal(names(jsonList$computations), c('computationDetails', 'pcoaVariance', 'computedVariable'))
  expect_equal(nrow(jsonList$computations), nMethods)
  expect_equal(jsonList$computations$pcoaVariance[[1]], c(15.3, 5.7))
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariable
  expect_equal(nrow(jsonList$computations$computedVariable), nMethods)
  expect_equal(names(jsonList$computations$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection','values'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$variableId[[1]], c('Axis1', 'Axis2'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$entityId[[1]], c('entity','entity'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataType[[1]], c('NUMBER', 'NUMBER'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataShape[[1]], c('CONTINUOUS', 'CONTINUOUS'))
  expect_equal(ncol(jsonList$computations$computedVariable$computedVariableDetails$values[[1]]), nrow(df))
  expect_false(jsonList$computations$computedVariable$computedVariableDetails$isCollection[[1]])
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariable$computedVariableMetadata), c('displayName'))
  expect_equal(jsonList$computations$computedVariable$computedVariableMetadata$displayName[[1]], c('Axis1 15.3%','Axis2 5.7%'))
  
  
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
  expect_equal(names(jsonList$computations), c('computationDetails', 'pcoaVariance', 'computedVariable'))
  expect_equal(nrow(jsonList$computations), nMethods)
  expect_equal(jsonList$computations$pcoaVariance[[1]], c(15.3, 5.7))
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariable
  expect_equal(nrow(jsonList$computations$computedVariable), nMethods)
  expect_equal(names(jsonList$computations$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection','values'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$variableId[[1]], c('Axis1', 'Axis2'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$entityId[[1]], c('entity','entity'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataType[[1]], c('NUMBER', 'NUMBER'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataShape[[1]], c('CONTINUOUS', 'CONTINUOUS'))
  expect_equal(ncol(jsonList$computations$computedVariable$computedVariableDetails$values[[1]]), nrow(df))
  expect_false(jsonList$computations$computedVariable$computedVariableDetails$isCollection[[1]])
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariable$computedVariableMetadata), c('displayName'))
  expect_equal(jsonList$computations$computedVariable$computedVariableMetadata$displayName[[1]], c('Axis1 15.3%','Axis2 5.7%'))

  
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
  
  dt <- betaDiv(df, "entity.SampleID", method='bray', verbose=T)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), 'entity.SampleID')
  expect_equal(typeof(dt$entity.SampleID), c('character'))
  attr <- attributes(dt)
  expect_equal(attr$computationDetails, "Error: beta diversity bray failed: missing values are not allowed with argument 'na.rm = FALSE'")
  expect_equal(typeof(attr$parameters), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$variableId), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$entityId), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$dataType), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$dataShape), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableMetadata$displayName), 'character')
  expect_equal(typeof(attr$pcoaVariance), 'double')
  
})