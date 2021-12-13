# Tests for ranked abundance functions
test_that('rankedAbundance returns a correctly formatted data.table', {
  
  df <- testOTU

  results <- rankedAbundance(df, "entity.SampleID", method='max', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), paste0('entity.',c('SampleID','Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
  results <- rankedAbundance(df, "entity.SampleID", method='median', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), paste0('entity.',c('SampleID','Lactobacillus','Snodgrassella','Gilliamella','Bifidobacterium','Frischella','Commensalibacter','unclassified Mitochondria','unclassified Rhizobiaceae','unclassified Chloroplast','Bombella')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
  results <- rankedAbundance(df, "entity.SampleID", method='variance', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), paste0('entity.',c('SampleID','Lactobacillus','Gilliamella','Bombella','Snodgrassella','Klebsiella','unclassified Rhizobiaceae','unclassified Enterobacterales','Fructobacillus','Pseudomonas','unclassified Chloroplast')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
  results <- rankedAbundance(df, "entity.SampleID", method='q3', verbose=F)
  expect_equal(nrow(results), nrow(df))
  expect_s3_class(results, 'data.table')
  expect_equal(names(results), paste0('entity.',c('SampleID','Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella')))
  expect_equal(unname(unlist(lapply(results, class))), c('character', rep('numeric',10)))
  
})


test_that("rankedAbundance returns a data.table with the correct attributes", {
  
  df <- testOTU
  
  results <- rankedAbundance(df, "entity.SampleID", method='median', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'median')
  expect_true(attr$isCutoff)
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('defaultRange','collectionVariable'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, c('Lactobacillus','Snodgrassella','Gilliamella','Bifidobacterium','Frischella','Commensalibacter','unclassified Mitochondria','unclassified Rhizobiaceae','unclassified Chloroplast','Bombella'))
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariable$computedVariableDetails$isCollection)
  expect_equal(attr$computedVariable$computedVariableMetadata$defaultRange, c(0, 1))
  
  results <- rankedAbundance(df, "entity.SampleID", method='max', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'max')
  expect_true(attr$isCutoff)
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('defaultRange','collectionVariable'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, c('Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella'))
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariable$computedVariableDetails$isCollection)
  expect_equal(attr$computedVariable$computedVariableMetadata$defaultRange, c(0, 1))
  
  results <- rankedAbundance(df, "entity.SampleID", method='q3', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'q3')
  expect_true(attr$isCutoff)
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('defaultRange','collectionVariable'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, c('Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella'))
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariable$computedVariableDetails$isCollection)
  expect_equal(attr$computedVariable$computedVariableMetadata$defaultRange, c(0, 1))
  
  results <- rankedAbundance(df, "entity.SampleID", method='variance', verbose=F)
  attr <- attributes(results)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'variance')
  expect_true(attr$isCutoff)
  expect_equal(names(attr$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
expect_equal(names(attr$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
  expect_equal(names(attr$computedVariable$computedVariableMetadata), c('defaultRange','collectionVariable'))
  expect_equal(attr$computedVariable$computedVariableDetails$variableId, c('Lactobacillus','Gilliamella','Bombella','Snodgrassella','Klebsiella','unclassified Rhizobiaceae','unclassified Enterobacterales','Fructobacillus','Pseudomonas','unclassified Chloroplast'))
  expect_equal(attr$computedVariable$computedVariableDetails$entityId, rep('entity',10))
  expect_true(attr$computedVariable$computedVariableDetails$isCollection)
  expect_equal(attr$computedVariable$computedVariableMetadata$defaultRange, c(0, 1))
  
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
  outJson <- getAppJson(appResults, "entity.SampleID")
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets','recordVariableDetails'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  # recordVariableDetails
  expect_equal(names(jsonList$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$recordVariableDetails$values), nrow(df))
  # computations
  expect_equal(names(jsonList$computations), c('computationDetails', 'isCutoff', 'computedVariable'))
  expect_equal(nrow(jsonList$computations), nMethods)
  expect_true(jsonList$computations$isCutoff[[1]])
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariable
  expect_equal(nrow(jsonList$computations$computedVariable), nMethods)
  expect_equal(names(jsonList$computations$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection','values'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$variableId[[1]], c('Lactobacillus','Snodgrassella','Gilliamella','Bifidobacterium','Frischella','Commensalibacter','unclassified Mitochondria','unclassified Rhizobiaceae','unclassified Chloroplast','Bombella'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$entityId[[1]], rep('entity', 10))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataType[[1]], rep('NUMBER', 10))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataShape[[1]], rep('CONTINUOUS', 10))
  expect_equal(ncol(jsonList$computations$computedVariable$computedVariableDetails$values[[1]]), nrow(df))
  expect_equal(nrow(jsonList$computations$computedVariable$computedVariableDetails$values[[1]]), 10)
  expect_true(jsonList$computations$computedVariable$computedVariableDetails$isCollection[[1]])
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariable$computedVariableMetadata), c('defaultRange','collectionVariable'))
  expect_equal(jsonList$computations$computedVariable$computedVariableMetadata$defaultRange[[1]], c(0, 1))
  expect_equal(jsonList$computations$computedVariable$computedVariableMetadata$collectionVariable$collectionType[[1]], 'abundance')
  
  # Supply only two methods
  nMethods <- 2
  appResults <- rankedAbundanceApp(df, "entity.SampleID", methods = c('q3', 'median'), verbose=F)
  outJson <- getAppJson(appResults, "entity.SampleID")
  jsonList <- jsonlite::fromJSON(outJson)
  expect_equal(names(jsonList), c('computations','parameterSets','recordVariableDetails'))
  expect_equal(length(jsonList$parameterSets), nMethods)
  # recordVariableDetails
  expect_equal(names(jsonList$recordVariableDetails), c('variableId','entityId','values'))
  expect_equal(jsonList$recordVariableDetails$variableId[[1]], 'SampleID')
  expect_equal(jsonList$recordVariableDetails$entityId[[1]], 'entity')
  expect_equal(length(jsonList$recordVariableDetails$values), nrow(df))
  # computations
  expect_equal(names(jsonList$computations), c('computationDetails', 'isCutoff', 'computedVariable'))
  expect_equal(nrow(jsonList$computations), nMethods)
  expect_true(jsonList$computations$isCutoff[[1]])
  # computationDetails
  expect_equal(length(jsonList$computations$computationDetails), nMethods)
  # computedVariable
  expect_equal(nrow(jsonList$computations$computedVariable), nMethods)
  expect_equal(names(jsonList$computations$computedVariable), c('computedVariableDetails','computedVariableMetadata'))
  # computedVariableDetails
  expect_equal(names(jsonList$computations$computedVariable$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection','values'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$variableId[[1]], c('Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella'))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$entityId[[1]], rep('entity', 10))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataType[[1]], rep('NUMBER', 10))
  expect_equal(jsonList$computations$computedVariable$computedVariableDetails$dataShape[[1]], rep('CONTINUOUS', 10))
  expect_equal(ncol(jsonList$computations$computedVariable$computedVariableDetails$values[[1]]), nrow(df))
  expect_equal(nrow(jsonList$computations$computedVariable$computedVariableDetails$values[[1]]), 10)
  expect_true(jsonList$computations$computedVariable$computedVariableDetails$isCollection[[1]])
  # computedVariableMetadata
  expect_equal(names(jsonList$computations$computedVariable$computedVariableMetadata), c('defaultRange','collectionVariable'))
  expect_equal(jsonList$computations$computedVariable$computedVariableMetadata$defaultRange[[1]], c(0, 1))
  expect_equal(jsonList$computations$computedVariable$computedVariableMetadata$collectionVariable$collectionType[[1]], 'abundance')
  
})

test_that("rankedAbundance results are consistent", {
  
  expect_snapshot_value({
    df <- testOTU
    appResults <- rankedAbundanceApp(df, "entity.SampleID", verbose=F)
  }, style = "serialize")
  
})