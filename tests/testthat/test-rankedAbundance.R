# Tests for ranked abundance functions
test_that('rankedAbundance returns a correctly formatted data.table', {
  
  df <- testOTU

  dt <- rankedAbundance(df, "entity.SampleID", method='max', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), paste0('entity.',c('SampleID','Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella')))
  expect_equal(unname(unlist(lapply(dt, class))), c('character', rep('numeric',10)))
  
  dt <- rankedAbundance(df, "entity.SampleID", method='median', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), paste0('entity.',c('SampleID','Lactobacillus','Snodgrassella','Gilliamella','Bifidobacterium','Frischella','Commensalibacter','unclassified Mitochondria','unclassified Rhizobiaceae','unclassified Chloroplast','Bombella')))
  expect_equal(unname(unlist(lapply(dt, class))), c('character', rep('numeric',10)))
  
  dt <- rankedAbundance(df, "entity.SampleID", method='variance', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), paste0('entity.',c('SampleID','Lactobacillus','Gilliamella','Bombella','Snodgrassella','Klebsiella','unclassified Rhizobiaceae','unclassified Enterobacterales','Fructobacillus','Pseudomonas','unclassified Chloroplast')))
  expect_equal(unname(unlist(lapply(dt, class))), c('character', rep('numeric',10)))
  
  dt <- rankedAbundance(df, "entity.SampleID", method='q3', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), paste0('entity.',c('SampleID','Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella')))
  expect_equal(unname(unlist(lapply(dt, class))), c('character', rep('numeric',10)))

  # With NAs
  nNAs <- 20
  df[sample(1:nrow(df), size=nNAs, replace = F), 2] <- NA
  expect_error(rankedAbundance(df, "entity.SampleID", method='q3', naToZero=F, verbose=F))
  dt <- rankedAbundance(df, "entity.SampleID", method='q3', verbose=F)
  expect_equal(nrow(dt), nrow(df))
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), paste0('entity.',c('SampleID','Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella')))
  expect_equal(unname(unlist(lapply(dt, class))), c('character', rep('numeric',10)))
  expect_equal(sum(dt > 0), 2847)   # ensure output is not all 0s.

  
})


test_that("rankedAbundance returns a data.table with the correct attributes", {
  
  df <- testOTU
  
  dt <- rankedAbundance(df, "entity.SampleID", method='median', verbose=F)
  attr <- attributes(dt)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'median')
  expect_true(attr$isCutoff)
  expect_equal(inherits(attr$computedVariable, "VariableMetadataList"), TRUE)
  expect_equal(unlist(lapply(as.list(attr$computedVariable[[1]]@members), function(x) {x@variableId})), c('Lactobacillus','Snodgrassella','Gilliamella','Bifidobacterium','Frischella','Commensalibacter','unclassified Mitochondria','unclassified Rhizobiaceae','unclassified Chloroplast','Bombella'))
  expect_equal(attr$computedVariable[[1]]@variableSpec@entityId, 'entity')
  expect_equal(attr$computedVariable[[1]]@displayRangeMin, 0)
  expect_equal(attr$computedVariable[[1]]@displayRangeMax, 1)
  
  dt <- rankedAbundance(df, "entity.SampleID", method='max', verbose=F)
  attr <- attributes(dt)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'max')
  expect_true(attr$isCutoff)
  expect_equal(inherits(attr$computedVariable, "VariableMetadataList"), TRUE)
  expect_equal(unlist(lapply(as.list(attr$computedVariable[[1]]@members), function(x) {x@variableId})), c('Gilliamella','Tyzzerella','Pseudomonas','Klebsiella','unclassified Chloroplast','Lactobacillus','unclassified Enterobacterales','Serratia','Frischella','Bombella'))
  expect_equal(attr$computedVariable[[1]]@variableSpec@entityId, 'entity')
  expect_equal(attr$computedVariable[[1]]@displayRangeMin, 0)
  expect_equal(attr$computedVariable[[1]]@displayRangeMax, 1)
  
  dt <- rankedAbundance(df, "entity.SampleID", method='q3', verbose=F)
  attr <- attributes(dt)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'q3')
  expect_true(attr$isCutoff)
  expect_equal(inherits(attr$computedVariable, "VariableMetadataList"), TRUE)
  expect_equal(unlist(lapply(as.list(attr$computedVariable[[1]]@members), function(x) {x@variableId})), c('Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella'))
  expect_equal(attr$computedVariable[[1]]@variableSpec@entityId, 'entity')
  expect_equal(attr$computedVariable[[1]]@displayRangeMin, 0)
  expect_equal(attr$computedVariable[[1]]@displayRangeMax, 1)
  
  dt <- rankedAbundance(df, "entity.SampleID", method='variance', verbose=F)
  attr <- attributes(dt)
  expect_true(all(c('computationDetails','parameters','computedVariable','isCutoff') %in% names(attr)))
  expect_equal(attr$parameters, 'variance')
  expect_true(attr$isCutoff)
  expect_equal(inherits(attr$computedVariable, "VariableMetadataList"), TRUE)
  expect_equal(unlist(lapply(as.list(attr$computedVariable[[1]]@members), function(x) {x@variableId})), c('Lactobacillus','Gilliamella','Bombella','Snodgrassella','Klebsiella','unclassified Rhizobiaceae','unclassified Enterobacterales','Fructobacillus','Pseudomonas','unclassified Chloroplast'))
  expect_equal(attr$computedVariable[[1]]@variableSpec@entityId, 'entity')
  expect_equal(attr$computedVariable[[1]]@displayRangeMin, 0)
  expect_equal(attr$computedVariable[[1]]@displayRangeMax, 1)
  
})

# no longer needed, tests exist in veupathUtils for writing VariableMetadata objects
#test_that("getMetadata() output is correctly represented in json", {
#  
#  df <- testOTU
#  
#  nMethods <- 1
#  results <- rankedAbundance(df, "entity.SampleID", method= 'q3', verbose=F)
#  outJson <- getMetadata(results)
#  jsonList <- jsonlite::fromJSON(outJson)
#  expect_equal(names(jsonList), c('computedVariableDetails','computedVariableMetadata'))
#  # computedVariableDetails
#  expect_equal(names(jsonList$computedVariableDetails), c('variableId','entityId','dataType','dataShape','isCollection'))
#  expect_equal(jsonList@variableSpec@variableId, c('Lactobacillus','Snodgrassella','Gilliamella','Frischella','Commensalibacter','unclassified Rhizobiaceae','Bifidobacterium','unclassified Mitochondria','unclassified Chloroplast','Bombella'))
#  expect_equal(jsonList@variableSpec@entityId, rep('entity', 10))
#  expect_equal(jsonList@variableSpec@dataType, rep('NUMBER', 10))
#  expect_equal(jsonList@variableSpec@dataShape, rep('CONTINUOUS', 10))
#  #expect_equal(ncol(jsonListi@variableSpec@values), nrow(df))
#  #expect_equal(nrow(jsonList@variableSpec@values), 10)
#  expect_true(jsonList@variableSpec@isCollection)
#  # computedVariableMetadata
#  expect_equal(names(jsonList$computedVariableMetadata), c('displayRangeMin','displayRangeMax','collectionVariable'))
#  expect_equal(jsonList@displayRangeMin, 0)
#  expect_equal(jsonList@displayRangeMax, 1)
#  
#})
#