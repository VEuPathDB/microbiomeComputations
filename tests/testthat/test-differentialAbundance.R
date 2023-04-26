# Tests for differential abundance methods
test_that('differentialAbundance returns a correctly formatted data.table', {
  df <- testOTU
  counts <- round(df[, -c("entity.SampleID")]*1000) # make into "counts"
  counts[ ,entity.SampleID:= df$entity.SampleID]
  nSamples <- dim(df)[1]
  sampleMetadata <- data.frame(list(
    "entity.SampleID" = df[["entity.SampleID"]],
    "entity.binA" = rep(c("binA_a", "binA_b"), nSamples/2, replace=T),
    "entity.cat2" = sample(c("cat2_a", "cat2_b"), nSamples, replace=T),
    "entity.cat3" = sample(paste0("cat3_", letters[1:3]), nSamples, replace=T),
    "entity.cat4" = sample(paste0("cat4_", letters[1:4]), nSamples, replace=T),
    "entity.contA" = rnorm(nSamples, sd=5)
    ))


  data <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts,
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')
  
  result <- differentialAbundance(data, comparisonVariable = "entity.binA", groupA = c('binA_a'), groupB = c('binA_b'), method='DESeq', verbose=F)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))


  ## Test with defined groups
  result <- differentialAbundance(data, comparisonVariable = "entity.cat3", groupA = c('cat3_a'), groupB = c('cat3_b', 'cat3_c'), method='DESeq', verbose=F)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))

  ## Expect diff abund to handle messy inputs

  # With no groups defined and a bin var

  # With samples ordered differently in the metadata and abundance data

})


test_that("differentialAbundance returns a data.table with the correct attributes" , {

  df <- testOTU
  counts <- round(df[, -c("entity.SampleID")]*1000) # make into "counts"
  counts[ ,entity.SampleID:= df$entity.SampleID]
  nSamples <- dim(df)[1]
  sampleMetadata <- data.frame(list(
    "entity.SampleID" = df[["entity.SampleID"]],
    "entity.binA" = sample(c("binA_a", "binA_b"), nSamples, replace=T),
    "entity.cat2" = sample(c("cat2_a", "cat2_b"), nSamples, replace=T),
    "entity.cat3" = sample(paste0("cat3_", letters[1:3]), nSamples, replace=T),
    "entity.cat4" = sample(paste0("cat4_", letters[1:4]), nSamples, replace=T)
    ))


  data <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts,
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')

  result <- differentialAbundance(data, comparisonVariable = "entity.binA", groupA = NULL, groupB = NULL, method='DESeq', verbose=F)
  expect_equal(result@parameters, 'comparisonVariable = entity.binA, groupA = binA_a, groupB = binA_b, method = DESeq')
  expect_equal(result@recordIdColumn, 'entity.SampleID')
})

test_that("differentialAbundance fails with improper inputs", {

  # Fail when we send in comparisonVar with > 2 values and no groups

  # Fail when group values are not found in the comparisonVar


})

test_that("differentialAbundance catches deseq errors", {

  df <- testOTU
  counts <- round(df[, -c("entity.SampleID")]*1000) # make into "counts"
  counts[ ,entity.SampleID:= df$entity.SampleID]
  nSamples <- dim(df)[1]
  sampleMetadata <- data.frame(list(
    "entity.SampleID" = df[["entity.SampleID"]],
    "entity.binA" = rep(c("binA_a", "binA_b"), nSamples/2, replace=T)
    ))

  # Use only a few taxa
  data <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts[, c("entity.SampleID","entity.1174-901-12","entity.A2")],
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')

  expect_error(differentialAbundance(data, comparisonVariable = "entity.binA", groupA = c('binA_a'), groupB = c('binA_b'), method='DESeq', verbose=T))


})