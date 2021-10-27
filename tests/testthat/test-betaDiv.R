# Tests for beta diversity functions
test_that('betaDiv returns something reasonable', {
  
  otu <- otu_local
  results <- betaDiv(otu, method='bray', verbose=F)
  expect_equal(dim(results$dt), c(288,133))
  
  
})