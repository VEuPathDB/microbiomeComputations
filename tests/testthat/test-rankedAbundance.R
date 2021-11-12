# Tests for ranked abundance functions
test_that('rankedAbundance returns a correctly formatted data.table', {
  
  df <- testOTU

  results <- rankedAbundance(df, "entity.SampleID", method='max', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_equal(names(results), paste0('entity.',c('SampleID','Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
  results <- rankedAbundance(df, "entity.SampleID", method='median', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_equal(names(results), paste0('entity.',c('SampleID','Lactobacillus','Snodgrassella','Gilliamella','Bifidobacterium','Frischella','Commensalibacter','unclassified Mitochondria','unclassified Rhizobiaceae','unclassified Chloroplast','Bombella')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
  results <- rankedAbundance(df, "entity.SampleID", method='variance', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_equal(names(results), paste0('entity.',c('SampleID','Lactobacillus','Gilliamella','Bombella','Snodgrassella','Klebsiella','unclassified Rhizobiaceae','unclassified Enterobacterales','Fructobacillus','Pseudomonas','unclassified Chloroplast')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
  results <- rankedAbundance(df, "entity.SampleID", method='q3', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_equal(names(results), paste0('entity.',c('SampleID','Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
})


test_that("rankedAbundance returns a data.table with the correct attributes", {
  
  df <- testOTU
  
  results <- rankedAbundance(df, "entity.SampleID", method='median', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails','isCutoff') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('variableId','entityId','dataType','dataShape','defaultRange', 'isCollection'))
  expect_equal(attr$computedVariableDetails$variableId, c('Lactobacillus','Snodgrassella','Gilliamella','Bifidobacterium','Frischella','Commensalibacter','unclassified Mitochondria','unclassified Rhizobiaceae','unclassified Chloroplast','Bombella'))
  expect_equal(attr$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariableDetails$isCollection)
  
  results <- rankedAbundance(df, "entity.SampleID", method='max', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails','isCutoff') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('variableId','entityId','dataType','dataShape','defaultRange', 'isCollection'))
  expect_equal(attr$computedVariableDetails$variableId, c('Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella'))
  expect_equal(attr$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariableDetails$isCollection)
  
  results <- rankedAbundance(df, "entity.SampleID", method='q3', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails','isCutoff') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('variableId','entityId','dataType','dataShape','defaultRange', 'isCollection'))
  expect_equal(attr$computedVariableDetails$variableId, c('Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella'))
  expect_equal(attr$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariableDetails$isCollection)
  
  results <- rankedAbundance(df, "entity.SampleID", method='variance', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','recordVariable','computedVariableDetails','isCutoff') %in% names(attr)))
  expect_equal(names(attr$computedVariableDetails), c('variableId','entityId','dataType','dataShape','defaultRange', 'isCollection'))
  expect_equal(attr$computedVariableDetails$variableId, c('Lactobacillus','Gilliamella','Bombella','Snodgrassella','Klebsiella','unclassified Rhizobiaceae','unclassified Enterobacterales','Fructobacillus','Pseudomonas','unclassified Chloroplast'))
  expect_equal(attr$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariableDetails$isCollection)
  
})


test_that("rankedAbundanceApp produces an appropriately structured list of computations", {
  
  df <- testOTU
  
  # Default - use all methods
  appResults <- rankedAbundanceApp(df, "entity.SampleID", verbose=F)
  expect_equal(length(appResults), 4)
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, ncol))), 11)
  expect_equal(unique(unlist(lapply(appResults, nrow))), nrow(df))
  
  # Test using a subset of methods
  appResults <- rankedAbundanceApp(df, "entity.SampleID", methods = c('q3', 'median'), verbose=F)
  expect_equal(length(appResults), 2)
  expect_equal(unique(unlist(lapply(appResults, class))), c('data.table','data.frame'))
  expect_equal(unique(unlist(lapply(appResults, ncol))), 11)
  expect_equal(unique(unlist(lapply(appResults, nrow))), nrow(df))
  
})

test_that("rankedAbundanceApp output is correctly represented in json", {
  
  df <- testOTU
  
  # Default all methods
  nMethods <- 4
  appResults <- rankedAbundanceApp(df, "entity.SampleID", verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','isCutoff', 'recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('variableId','entityId','dataType','dataShape','defaultRange','isCollection','values'))
  expect_equal(nrow(jsonList$computations$computedVariableDetails), nMethods)
  expect_equal(jsonList$computations$computedVariableDetails$variableId[[2]], c('Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella'))
  expect_equal(jsonList$computations$computedVariableDetails$entityId[[1]], rep('entity',10))
  expect_equal(jsonList$computations$computedVariableDetails$defaultRange[[1]], c(0, 1))
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  expect_equal(nrow(jsonList$computations$recordVariableDetails), nMethods)
  expect_equal(names(jsonList$computations$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$computations$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$computations$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$computations$recordVariableDetails$values[[1]]), nrow(df))
  expect_equal(length(jsonList$computations$isCutoff), nMethods)
  expect_equal(jsonList$computations$isCutoff, rep(TRUE, nMethods))
  
  # Supply only two methods
  nMethods <- 2
  appResults <- rankedAbundanceApp(df, "entity.SampleID", methods = c('q3', 'median'), verbose=F)
  outJson <- getAppJson(appResults)
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  expect_equal(names(jsonList$computations), c('computedVariableDetails','computationDetails','isCutoff', 'recordVariableDetails'))
  expect_equal(names(jsonList$computations$computedVariableDetails), c('variableId','entityId','dataType','dataShape','defaultRange','isCollection','values'))
  expect_equal(nrow(jsonList$computations$computedVariableDetails), nMethods)
  expect_equal(jsonList$computations$computedVariableDetails$variableId[[1]], c('Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella'))
  expect_equal(jsonList$computations$computedVariableDetails$entityId[[1]], rep('entity',10))
  expect_equal(jsonList$computations$computedVariableDetails$defaultRange[[1]], c(0, 1))
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  expect_equal(nrow(jsonList$computations$recordVariableDetails), nMethods)
  expect_equal(names(jsonList$computations$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$computations$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$computations$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$computations$recordVariableDetails$values[[1]]), nrow(df))
  expect_equal(length(jsonList$computations$isCutoff), nMethods)
  expect_equal(jsonList$computations$isCutoff, rep(TRUE, nMethods))
})

test_that("rankedAbundance results are consistent", {
  
  expect_snapshot_value({
    df <- testOTU
    appResults <- rankedAbundanceApp(df, "entity.SampleID", verbose=F)
  }, style = "serialize")
  
})