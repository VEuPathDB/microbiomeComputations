# Tests for alpha diversity functions
test_that('alphaDiv returns something reasonable', {
  
  otu <- otu_local
  results <- alphaDiv(otu, method='shannon', verbose=F)
  expect_equal(dim(results$dt), c(288,2))
  
  
})