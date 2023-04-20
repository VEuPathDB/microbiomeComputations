test_that('ComputeResult validation works', {
  
  df <- testOTU

  expect_error(microbiomeComputations::AbundanceData(
              data = df))
  expect_error(microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = c('entity.SampleID', 'test')))
  expect_error(microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = c('test')))
  expect_error(microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = c('entity.SampleID'),
              ancestorIdColumns = c('test')))

  df$entity.strings <- 'a'

  expect_error(microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = c('entity.SampleID')))

  df$entity.strings <- .1
  names(df)[names(df) == 'entity.strings'] <- 'entity2.notstrings'

  expect_error(microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = c('entity.SampleID')))

  df$entity2.notstrings <- NULL
  testing <- microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = c('entity.SampleID'))

  expect_true(inherits(testing, 'AbundanceData'))
  expect_equal(slotNames(testing), c('data', 'sampleMetadata', 'recordIdColumn', 'ancestorIdColumns', 'imputeZero'))
  expect_equal(nrow(testing@data), 288)
  expect_equal(ncol(testing@data), 909)
  expect_equal(testing@recordIdColumn, 'entity.SampleID')
  expect_equal(testing@ancestorIdColumns, character(0))
  expect_true(testing@imputeZero)
})