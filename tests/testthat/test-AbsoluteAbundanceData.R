test_that('AbsoluteAbundanceData validation works', {

  # Most tests are handled by the parent class.

  df <- testOTU

  # Expect error when input data is not integers (not rounded)
  expect_error(microbiomeComputations::AbsoluteAbundanceData(
              data = df,
              recordIdColumn = c('entity.SampleID')))


  counts <- round(df[, -c("entity.SampleID")]*1000) # make into "counts"
  counts[ ,entity.SampleID:= df$entity.SampleID]
  testing <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts,
              recordIdColumn = c('entity.SampleID'))

  expect_true(inherits(testing, 'AbundanceData'))
  expect_true(inherits(testing, 'AbsoluteAbundanceData'))
  expect_equal(slotNames(testing), c('data', 'sampleMetadata', 'recordIdColumn', 'ancestorIdColumns', 'imputeZero', 'removeEmptySamples'))
  expect_equal(nrow(testing@data), 288)
  expect_equal(ncol(testing@data), 909)
  expect_equal(testing@recordIdColumn, 'entity.SampleID')
  expect_equal(testing@ancestorIdColumns, character(0))
  expect_true(testing@imputeZero)
})