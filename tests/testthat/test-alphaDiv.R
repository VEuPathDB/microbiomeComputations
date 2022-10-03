# Tests for alpha diversity functions
test_that('alphaDiv returns a correctly formatted data.table', {

  df <- testOTU
  
  dt <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(class(dt$entity.SampleID), class(dt$entity.alphaDiversity)), c('character','numeric'))
  
  dt <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(class(dt$entity.SampleID), class(dt$entity.alphaDiversity)), c('character','numeric'))
  
  dt <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(class(dt$entity.SampleID), class(dt$entity.alphaDiversity)), c('character','numeric'))

  # With NAs
  nNAs <- 20
  df[sample(1:nrow(df), size=nNAs, replace = F), 2] <- NA
  dt <- alphaDiv(df, "entity.SampleID", method='shannon', naToZero = F, verbose=F)  # vegan diversity sets output to NA if input has NA. See issue #187
  expect_equal(sum(is.na(dt)), nNAs)
  dt <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  expect_equal(sum(is.na(dt)), 0)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID','entity.alphaDiversity'))
  expect_equal(c(class(dt$entity.SampleID), class(dt$entity.alphaDiversity)), c('character','numeric'))
  expect_equal(sum(dt>0), 576)  # ensure output is not all 0s.
})

test_that("alphaDiv returns a data.table with the correct attributes", {
  
  df <- testOTU
  
  dt <- alphaDiv(df, "entity.SampleID", method='shannon', verbose=F)
  attr <- attributes(dt)
  expect_true(all(c('computationDetails','parameters','computedVariable') %in% names(attr)))
  expect_equal(attr$parameters, 'shannon')
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('displayName','displayRangeMin','displayRangeMax'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayName, 'Shannon Diversity')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayRangeMin, '0')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayRangeMax, '1')
  
  dt <- alphaDiv(df, "entity.SampleID", method='simpson', verbose=F)
  attr <- attributes(dt)
  expect_true(all(c('computationDetails','parameters','computedVariable') %in% names(attr)))
  expect_equal(attr$parameters, 'simpson')
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('displayName','displayRangeMin','displayRangeMax'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayName, 'Simpson Diversity')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayRangeMin, '0')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayRangeMax, '1')
  
  dt <- alphaDiv(df, "entity.SampleID", method='evenness', verbose=F)
  attr <- attributes(dt)
  expect_true(all(c('computationDetails','parameters','computedVariable') %in% names(attr)))
  expect_equal(attr$parameters, 'evenness')
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('displayName','displayRangeMin','displayRangeMax'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, 'entity')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayName, 'Pielou\'s Evenness')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayRangeMin, '0')
  expect_equal(attr$computedVariable$computedVariableMetadata$displayRangeMax, '1')
  
})


test_that("getMetadata() output is correctly represented in json", {
  
  df <- testOTU
  
  nMethods <- 1
  results <- alphaDiv(df, 'entity.SampleID', method='evenness', verbose=F)
  outJson <- getMetadata(results)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computedVariableDetails), c('variableId','entityId','dataType','dataShape'))
  expect_equal(jsonList$computedVariableDetails$variableId, 'alphaDiversity')
  expect_equal(jsonList$computedVariableDetails$entityId, 'entity')
  expect_equal(jsonList$computedVariableDetails$dataType, 'NUMBER')
  expect_equal(jsonList$computedVariableDetails$dataShape, 'CONTINUOUS')
  #expect_equal(length(jsonList$computedVariableDetails$values), nrow(df))
  # computedVariableMetadata
  expect_equal(names(jsonList$computedVariableMetadata), c('displayName','displayRangeMin','displayRangeMax'))
  expect_equal(jsonList$computedVariableMetadata$displayName, 'Pielou\'s Evenness')
  expect_equal(jsonList$computedVariableMetadata$displayRangeMin, '0')
  expect_equal(jsonList$computedVariableMetadata$displayRangeMax, '1')
})


test_that("alphaDiv fails gracefully", {
  
  df <- testOTU
  df$entity.Abiotrophia <- NA
  
  dt <- alphaDiv(df, "entity.SampleID", method='simpson', naToZero=F, verbose=T)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('entity.SampleID', 'entity.alphaDiversity'))
  expect_true(all(is.na(dt[['entity.alphaDiversity']])))
  expect_equal(typeof(dt$entity.SampleID), c('character'))
  attr <- attributes(dt)
  # expect_equal(attr$computationDetails, "Error: alpha diversity evenness failed: input data must be non-negative")
  expect_equal(typeof(attr$parameters), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$variableId), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$entityId), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$dataType), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableDetails$dataShape), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableMetadata$displayName), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableMetadata$displayRangeMin), 'character')
  expect_equal(typeof(attr$computedVariable$computedVariableMetadata$displayRangeMax), 'character')
  
  
})
