# Tests for beta diversity functions
test_that('betaDiv returns something reasonable', {
  
  otu <- testOTU
  results <- betaDiv(otu, "SampleID", method='bray', verbose=F)
  expect_equal(NROW(results), 288)
  
  results <- betaDiv(otu, "SampleID", method='jaccard', verbose=F)
  expect_equal(NROW(results), 288)
  
  results <- betaDiv(otu, "SampleID", method='jsd', verbose=F)
  expect_equal(NROW(results), 288)
  
  
})

test_that("betaDivApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- betaDivApp(otu, "SampleID", methods=c('bray'), verbose=F)
  expect_equal(length(appResults), 1)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 1)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','pcoaVariance','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','isCollection','values'))
  
  appResults <- betaDivApp(otu, "SampleID", methods = c('bray','jaccard'), verbose=F)
  expect_equal(length(appResults), 2)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 2)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','pcoaVariance','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','isCollection','values'))
  
})