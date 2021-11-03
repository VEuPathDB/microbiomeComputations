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
  
  appResults <- alphaDivApp(otu, F)
  expect_equal(names(appResults), c('shannon','simpson','evenness'))
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('shannon','simpson','evenness'))
  expect_equal(names(jsonList$shannon), c('data','computedVariableDetails','parameterSet','computationDetails'))
  expect_equal(names(jsonList$shannon$computedVariableDetails), c('id','entity','displayLabel','defaultRange'))
})