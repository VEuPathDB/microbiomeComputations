# Tests for ranked abundance functions
test_that('rankedAbundance returns something reasonable', {
  
  otu <- testOTU
  results <- rankedAbundance(otu, method='max', verbose=F)
  expect_equal(NROW(results), 288)
  
  
})


test_that("rankedAbundanceApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- rankedAbundanceApp(otu, verbose=F)
  expect_equal(length(appResults), 4)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('data','computedVariableDetails','parameters','computationDetails','isCutoff'))
  expect_equal(names(jsonList$computedVariableDetails), c('id','entity','defaultRange','isCollection'))
  
  appResults <- rankedAbundanceApp(otu, methods = c('q3', 'median'), verbose=F)
  expect_equal(length(appResults), 2)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('data','computedVariableDetails','parameters','computationDetails','isCutoff'))
  expect_equal(names(jsonList$computedVariableDetails), c('id','entity','defaultRange','isCollection'))
  
})