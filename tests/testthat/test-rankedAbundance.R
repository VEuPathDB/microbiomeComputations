# Tests for ranked abundance functions
test_that('rankedAbundance returns something reasonable', {
  
  otu <- testOTU
  results <- rankedAbundance(otu, "SampleID", method='max', verbose=F)
  expect_equal(NROW(results), 288)
  expect_equal(names(results), c('SampleID','Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella'))
  
  
})


test_that("rankedAbundanceApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- rankedAbundanceApp(otu, "SampleID", verbose=F)
  expect_equal(length(appResults), 4)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 4)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','isCutoff','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','defaultRange','isCollection','values'))
  
  appResults <- rankedAbundanceApp(otu, "SampleID", methods = c('q3', 'median'), verbose=F)
  expect_equal(length(appResults), 2)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 2)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','isCutoff','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','defaultRange','isCollection','values'))
  
})