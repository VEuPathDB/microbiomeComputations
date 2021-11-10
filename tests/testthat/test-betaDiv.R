# Tests for beta diversity functions
test_that('betaDiv returns a correctly formatted data.table', {
  
  df <- testOTU

  results <- betaDiv(df, "entity.SampleID", method='bray', verbose=F)
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(results, class))), c('character','numeric','numeric'))
  
  results <- betaDiv(df, "entity.SampleID", method='jaccard', verbose=F)
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(results, class))), c('character','numeric','numeric'))
  
  results <- betaDiv(df, "entity.SampleID", method='jsd', verbose=F)
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unname(unlist(lapply(results, class))), c('character','numeric','numeric'))
  
  
})

test_that("betaDiv returns a data.table with the correct attributes" , {
  
  df <- testOTU
  
  results <- betaDiv(df, "entity.SampleID", method='bray', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails','pcoaVariance') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('id','entity','displayLabel','isCollection'))
  expect_equal(attr$computedVariableDetails$id, c('Axis1','Axis2'))
  expect_equal(attr$computedVariableDetails$entity, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, c('Axis1 15.3%','Axis2 5.7%'))
  
  results <- betaDiv(df, "entity.SampleID", method='jaccard', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails','pcoaVariance') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('id','entity','displayLabel','isCollection'))
  expect_equal(attr$computedVariableDetails$id, c('Axis1','Axis2'))
  expect_equal(attr$computedVariableDetails$entity, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, c('Axis1 10.0%','Axis2 4.3%'))
  
  results <- betaDiv(df, "entity.SampleID", method='jsd', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails','pcoaVariance') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('id','entity','displayLabel','isCollection'))
  expect_equal(attr$computedVariableDetails$id, c('Axis1','Axis2'))
  expect_equal(attr$computedVariableDetails$entity, 'entity')
  expect_equal(attr$computedVariableDetails$displayLabel, c('Axis1 25.2%','Axis2 17.5%'))
})


test_that("betaDivApp produces an appropriately structured list of computations", {
  
  df <- testOTU
  
  # Default - use all methods
  appResults <- betaDivApp(df, "entity.SampleID", verbose=F)
  expect_equal(length(appResults), 3)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, NROW))), NROW(df))
  
  # Test using a subset of methods
  appResults <- betaDivApp(df, "entity.SampleID", methods = c('jaccard','bray'), verbose=F)
  expect_equal(length(appResults), 2)
  expect_equal(unique(unlist(lapply(appResults,names))), c('entity.SampleID','entity.Axis1','entity.Axis2'))
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, NROW))), NROW(df))
  
})


test_that("betaDivApp output is correctly represented in json", {
  
  df <- testOTU
  
  appResults <- betaDivApp(df, "entity.SampleID", methods=c('bray'), verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 1)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','pcoaVariance','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','isCollection','values'))
  expect_equal(jsonList$computations$computedVariableDetails$id[[1]], c('Axis1','Axis2'))
  
  appResults <- betaDivApp(df, "entity.SampleID", methods = c('bray','jaccard'), verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 2)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','pcoaVariance','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','isCollection','values'))
  expect_equal(jsonList$computations$computedVariableDetails$id[[1]], c('Axis1','Axis2'))
  
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
  expect_equal(NROW(results), NROW(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), 'entity.SampleID')
  expect_equal(typeof(results$entity.SampleID), c('character'))
  attr <- attributes(results)
  expect_equal(attr$computationDetails, "Error: beta diversity bray failed: missing values are not allowed with argument 'na.rm = FALSE'")
  expect_equal(typeof(attr$parameters), 'character')
  expect_equal(typeof(attr$recordVariable), 'character')
  expect_equal(typeof(attr$pcoaVariance), 'double')
  expect_equal(typeof(attr$computedVariableDetails$id), 'character')
  expect_equal(typeof(attr$computedVariableDetails$entity), 'character')
  expect_equal(typeof(attr$computedVariableDetails$displayLabel), 'character')
  expect_equal(typeof(attr$computedVariableDetails$isCollection), 'logical')
  
  
  
})