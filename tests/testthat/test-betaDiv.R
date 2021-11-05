# Tests for beta diversity functions
test_that('betaDiv returns something reasonable', {
  
  otu <- testOTU
  results <- betaDiv(otu, "entity.SampleID", method='bray', verbose=F)
  expect_equal(NROW(results), 288)
  expect_equal(names(results), c('entity.SampleID','entity.Axis.1','entity.Axis.2'))
  
  results <- betaDiv(otu, "entity.SampleID", method='jaccard', verbose=F)
  expect_equal(NROW(results), 288)
  expect_equal(names(results), c('entity.SampleID','entity.Axis.1','entity.Axis.2'))
  
  results <- betaDiv(otu, "entity.SampleID", method='jsd', verbose=F)
  expect_equal(NROW(results), 288)
  expect_equal(names(results), c('entity.SampleID','entity.Axis.1','entity.Axis.2'))
  
  
})

test_that("betaDivApp doesn't fail", {
  
  otu <- testOTU
  
  appResults <- betaDivApp(otu, "entity.SampleID", methods=c('bray'), verbose=F)
  expect_equal(length(appResults), 1)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 1)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','pcoaVariance','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','isCollection','values'))
  expect_equal(jsonList$computations$computedVariableDetails$id[[1]], c('Axis.1','Axis.2'))
  
  appResults <- betaDivApp(otu, "entity.SampleID", methods = c('bray','jaccard'), verbose=F)
  expect_equal(length(appResults), 2)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(NROW(jsonList$computations), 2)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','pcoaVariance','recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('id','entity','displayLabel','isCollection','values'))
  expect_equal(jsonList$computations$computedVariableDetails$id[[1]], c('Axis.1','Axis.2'))
  
})