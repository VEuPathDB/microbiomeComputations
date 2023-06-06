# Tests for differential abundance methods

test_that('differentialAbundance returns a correctly formatted data.table', {
  df <- data.table::copy(testOTU)
  counts_test <- round(testOTU[, -c("entity.SampleID")]*1000) # make into "counts"
  counts_test[ ,entity.SampleID:= testOTU$entity.SampleID]
  nSamples <- dim(testOTU)[1]
  sampleMetadata_test <- data.frame(list(
    "entity.SampleID" = testOTU[["entity.SampleID"]],
    "entity.binA" = rep(c("binA_a", "binA_b"), nSamples/2, replace=T),
    "entity.cat3" = sample(paste0("cat3_", letters[1:3]), nSamples, replace=T),
    "entity.cat4" = sample(paste0("cat4_", letters[1:4]), nSamples, replace=T)
    ))


  data <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts_test,
              sampleMetadata = sampleMetadata_test,
              recordIdColumn = 'entity.SampleID')
  
  # Binary comparisonVariable
  result <- differentialAbundance(data, comparisonVariable = "entity.binA", groupA = c('binA_a'), groupB = c('binA_b'), method='DESeq', verbose=F)
  expect_equal(length(result@droppedColumns), 182)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))
  expect_true(all(!is.na(stats[, c('log2foldChange', 'pValue', 'pointID')])))


  ## Test with defined groups
  result <- differentialAbundance(data, comparisonVariable = "entity.cat3", groupA = c('cat3_a'), groupB = c('cat3_b', 'cat3_c'), method='DESeq', verbose=F)
  expect_equal(length(result@droppedColumns), 182)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))
  expect_true(all(!is.na(stats[, c('log2foldChange', 'pValue', 'pointID')])))


  # When defined groups end up subsetting the incoming data
  result <- differentialAbundance(data, comparisonVariable = "entity.cat4", groupA = c('cat4_a'), groupB = c('cat4_b'), method='DESeq', verbose=F)
  expect_equal(length(result@droppedColumns), 381)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  expect_equal(nrow(sampleMetadata_test[entity.cat4 %in% c('cat4_a','cat4_b'), ]), nrow(dt))
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))
  expect_true(all(!is.na(stats[, c('log2foldChange', 'pValue', 'pointID')])))




  ## Expect diff abund to handle messy inputs

  # With no groups defined and a bin var
  result <- differentialAbundance(data, comparisonVariable = "entity.binA", groupA = NULL, groupB = NULL, method='DESeq', verbose=F)
  expect_equal(length(result@droppedColumns), 182)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))
  expect_true(all(!is.na(stats[, c('log2foldChange', 'pValue', 'pointID')])))


  # With only some comparisonVariable values found in the metadata
  result <- differentialAbundance(data, comparisonVariable = "entity.cat4", groupA = c('cat4_a','test'), groupB = c('cat4_b'), method='DESeq', verbose=F)
  expect_equal(length(result@droppedColumns), 381)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  expect_equal(nrow(sampleMetadata_test[entity.cat4 %in% c('cat4_a','cat4_b'), ]), nrow(dt))
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))
  expect_true(all(!is.na(stats[, c('log2foldChange', 'pValue', 'pointID')])))

  # With samples ordered differently in the metadata and abundance data
  data@data <- data@data[288:1, ]
  result <- differentialAbundance(data, comparisonVariable = "entity.cat4", groupA = c('cat4_a'), groupB = c('cat4_b'), method='DESeq', verbose=T)
  expect_equal(length(result@droppedColumns), 381)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))
  expect_true(all(!is.na(stats[, c('log2foldChange', 'pValue', 'pointID')])))

})


test_that("differentialAbundance returns a ComputeResult with the correct slots" , {

  df <- data.table::copy(testOTU)
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
  expect_equal(class(result@droppedColumns), 'character')
})

test_that("differentialAbundance fails with improper inputs", {

  df <- data.table::copy(testOTU)
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

  # Fail when we send in comparisonVar with > 2 values and do not provide one of groupA/B
  expect_error(differentialAbundance(data, comparisonVariable = "entity.cat4", groupA = NULL, groupB = c('cat4a','cat4c'), method='DESeq', verbose=F))

  # Fail when group values are not found in the comparisonVar
  expect_error(differentialAbundance(data, comparisonVariable = "entity.cat4", groupA = c('c','cat4_b'), groupB = c('a','b'), method='DESeq', verbose=F))

  # Fail when values are duplicated in groupA and groupB
  expect_error(differentialAbundance(data, comparisonVariable = "entity.cat4", groupA = c('cat4_a','cat4_b'), groupB = c('cat4_a','cat4_c'), method='DESeq', verbose=F))

})

test_that("differentialAbundance catches deseq errors", {

  df <- data.table::copy(testOTU)
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