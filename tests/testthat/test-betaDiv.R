# Tests for beta diversity functions
test_that('betaDiv returns something reasonable', {
  
  otu <- testOTU
  results <- betaDiv(otu, method='bray', verbose=F)
  expect_equal(dim(results$dt), c(288,133))
  
  results <- betaDiv(otu, method='jsd', verbose=T)
  
  
})