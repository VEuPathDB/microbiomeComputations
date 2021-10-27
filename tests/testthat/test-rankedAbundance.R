# Tests for ranked abundance functions
test_that('rankedAbundance returns something reasonable', {
  
  otu <- otu_local
  results <- rankedAbundance(otu, method='max', verbose=F)
  expect_equal(length(results$computedVariables), 10)
  
  
})