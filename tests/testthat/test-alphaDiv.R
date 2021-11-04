# Tests for alpha diversity functions
test_that('alphaDiv returns something reasonable', {
  
  otu <- testOTU
  results <- alphaDiv(otu, method='shannon', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('record','alphaDiversity'))
  expect_equal(c(typeof(results$record), typeof(results$alphaDiversity)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Shannon')
  
  results <- alphaDiv(otu, method='simpson', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('record','alphaDiversity'))
  expect_equal(c(typeof(results$record), typeof(results$alphaDiversity)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Simpson')
  
  results <- alphaDiv(otu, method='evenness', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('record','alphaDiversity'))
  expect_equal(c(typeof(results$record), typeof(results$alphaDiversity)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariableDetails$id, 'alphaDiversity')
  expect_equal(attr$computedVariableDetails$displayLabel, 'Pielou\'s Evenness')
  
})

test_that("alphaDivApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- alphaDivApp(otu, verbose=F)
  expect_equal(length(appResults), 3)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(NROW(jsonList), 3)
  expect_equal(names(jsonList), c('data','computedVariableDetails','parameters','computationDetails'))
  expect_equal(names(jsonList$computedVariableDetails), c('id','entity','displayLabel','defaultRange'))
  
  appResults <- alphaDivApp(otu, methods=c('evenness','simpson'), verbose=F)
  expect_equal(length(appResults), 2)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(NROW(jsonList), 2)
  expect_equal(names(jsonList), c('data','computedVariableDetails','parameters','computationDetails'))
  expect_equal(names(jsonList$computedVariableDetails), c('id','entity','displayLabel','defaultRange'))
})