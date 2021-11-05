# Tests for alpha diversity functions
test_that('alphaDiv returns something reasonable', {
  
  otu <- testOTU
  results <- alphaDiv(otu, "SampleID", method='shannon', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('SampleID','alphaDiversity'))
  expect_equal(c(typeof(results$SampleID), typeof(results$alphaDiversity)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Shannon')
  
  results <- alphaDiv(otu, "SampleID", method='simpson', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('SampleID','alphaDiversity'))
  expect_equal(c(typeof(results$SampleID), typeof(results$alphaDiversity)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Simpson')
  
  results <- alphaDiv(otu, "SampleID", method='evenness', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('SampleID','alphaDiversity'))
  expect_equal(c(typeof(results$SampleID), typeof(results$alphaDiversity)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Pielou\'s Evenness')
  
})

test_that("alphaDivApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- alphaDivApp(otu, "SampleID", verbose=F)
  expect_equal(length(appResults), 3)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 3)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','defaultRange','values'))
  
  appResults <- alphaDivApp(otu, "SampleID", methods=c('evenness','simpson'), verbose=F)
  expect_equal(length(appResults), 2)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 2)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','defaultRange','values'))
})