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


# after refactoring to a parent 'computation' class this should move somewhere else
test_that("getMetadata() output is correctly represented in json", {
  
  df <- testOTU
  
  nMethods <- 1
  results <- betaDiv(df, "entity.SampleID", method=c('bray'), verbose=F)
  outJson <- getMetadata(results)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(jsonList$computedVariableDetails$variableId, c('Axis1', 'Axis2'))
  expect_equal(jsonList$computedVariableDetails$entityId, c('entity','entity'))
  expect_equal(jsonList$computedVariableDetails$dataType, c('NUMBER', 'NUMBER'))
  expect_equal(jsonList$computedVariableDetails$dataShape, c('CONTINUOUS', 'CONTINUOUS'))
  expect_false(jsonList$computedVariableDetails$isCollection)
  # computedVariableMetadata
  expect_equal(names(jsonList$computedVariableMetadata), c('displayName'))
  expect_equal(jsonList$computedVariableMetadata$displayName, c('Axis1 15.3%','Axis2 5.7%'))

})


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
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$variableId), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$entityId), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$dataType), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$dataShape), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableMetadata$displayName), 'character')
  expect_equal(typeof(attr$pcoaVariance), 'double')
  
})