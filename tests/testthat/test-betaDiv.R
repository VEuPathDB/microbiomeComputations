# Tests for beta diversity functions
test_that('betaDiv returns something reasonable', {
  
  otu <- testOTU
  results <- betaDiv(otu, method='bray', verbose=F)
  expect_equal(NROW(results), 288)
  
  results <- betaDiv(otu, method='jaccard', verbose=F)
  expect_equal(NROW(results), 288)
  
  results <- betaDiv(otu, method='jsd', verbose=F)
  expect_equal(NROW(results), 288)
  
  
})

test_that("betaDivApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- betaDivApp(otu, 'bray', verbose=F)
  expect_equal(names(appResults), c('bray'))
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('bray'))
  expect_equal(names(jsonList$bray), c('data','computedVariableDetails','parameterSet','computationDetails','pcoaVariance'))
  
})