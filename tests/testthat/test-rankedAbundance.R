# Tests for ranked abundance functions
test_that('rankedAbundance returns something reasonable', {
  
  otu <- testOTU
  results <- rankedAbundance(otu, method='max', verbose=F)
  expect_equal(NROW(results), 288)
  
  
})


test_that("rankedAbundanceApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- rankedAbundanceApp(otu, verbose=F)
  expect_equal(names(appResults), c('median','max','q3','var'))
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('median','max','q3','var'))
  expect_equal(names(jsonList$median), c('data','computedVariableDetails','parameterSet','computationDetails','isCutoff'))
  
})