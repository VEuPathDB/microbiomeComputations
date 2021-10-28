# Tests for alpha diversity functions
test_that('alphaDiv returns something reasonable', {
  
  otu <- testOTU
  results <- alphaDiv(otu, method='shannon', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('SampleID','shannon'))
  expect_equal(c(typeof(results$SampleID), typeof(results$shannon)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariables, 'shannon')
  expect_equal(attr$computedVariableLabels, 'Shannon')
  
  results <- alphaDiv(otu, method='simpson', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('SampleID','simpson'))
  expect_equal(c(typeof(results$SampleID), typeof(results$simpson)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariables, 'simpson')
  expect_equal(attr$computedVariableLabels, 'Simpson')
  
  results <- alphaDiv(otu, method='evenness', verbose=F)
  expect_equal(NROW(results), 288)
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), c('SampleID','evenness'))
  expect_equal(c(typeof(results$SampleID), typeof(results$evenness)), c('character','double'))
  attr <- attributes(results)
  expect_equal(attr$computedVariables, 'evenness')
  expect_equal(attr$computedVariableLabels, 'Pielou\'s Evenness')
  
})

test_that("alphaDivApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- alphaDivApp(otu, F)
  expect_equal(names(appResults), c('shannon','simpson','evenness'))
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('shannon','simpson','evenness'))
  expect_equal(names(jsonList$shannon), c('data','metaData'))
  expect_equal(names(jsonList$shannon$metaData), c('computeDetails','defaultRange','computedAxisLabel','computedVariables','computedVariableLabels'))
})