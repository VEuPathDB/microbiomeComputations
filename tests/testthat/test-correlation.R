# Tests for differential abundance methods

test_that('correlation returns a correctly formatted data.table', {
  df <- testOTU
  nSamples <- dim(df)[1]
  sampleMetadata <- data.frame(list(
    "entity.SampleID" = df[["entity.SampleID"]],
    "entity.contA" = rnorm(nSamples),
    "entity.contB" = rnorm(nSamples),
    "entity.contC" = rnorm(nSamples),
    "entity.dateA" = sample(seq(as.Date('1999/01/01'), as.Date('2000/01/01'), by="day"), nSamples)
    ))


  data <- microbiomeComputations::AbundanceData(
              data = df,
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')
  
  ## All numeric sample variables
  result <- correlation(data, 'pearson', variables = NULL, verbose = FALSE)
  # Check data (only sample ids)
  dt <- result@data
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('SampleID'))
  expect_equal(nrow(dt), nSamples)
  # Check stats (all correlation outputs)
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('var1','var2','correlationCoef'))
  expect_equal(nrow(stats), (ncol(testOTU) - 1) * length(veupathUtils::findNumericCols(sampleMetadata))) # Should be number of taxa * number of metadata vars
  expect_equal(as.character(unique(stats$var1)), names(testOTU)[2:length(names(testOTU))])
  expect_equal(as.character(unique(stats$var2)), c('entity.contA', 'entity.contB', 'entity.contC'))


  ## With method = spearman
  result <- correlation(data, 'spearman', variables = NULL, verbose = FALSE)
  # Check data (only sample ids)
  dt <- result@data
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('SampleID'))
  expect_equal(nrow(dt), nSamples)
  # Check stats (all correlation outputs)
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('var1','var2','correlationCoef'))
  expect_equal(nrow(stats), (ncol(testOTU) - 1) * length(veupathUtils::findNumericCols(sampleMetadata))) # Should be number of taxa * number of metadata vars
  expect_equal(as.character(unique(stats$var1)), names(testOTU)[2:length(names(testOTU))])
  expect_equal(as.character(unique(stats$var2)), c('entity.contA', 'entity.contB', 'entity.contC'))


  ## With samples ordered differently in the abundance data and metadata
  data@data <- data@data[288:1, ]
  result <- correlation(data, 'spearman', variables = NULL, verbose = FALSE)
  # Check data (only sample ids)
  dt <- result@data
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('SampleID'))
  expect_equal(nrow(dt), nSamples)
  # Check stats (all correlation outputs)
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('var1','var2','correlationCoef'))
  expect_equal(nrow(stats), (ncol(testOTU) - 1) * length(veupathUtils::findNumericCols(sampleMetadata))) # Should be number of taxa * number of metadata vars
  expect_equal(as.character(unique(stats$var1)), names(testOTU)[2:length(names(testOTU))])
  expect_equal(as.character(unique(stats$var2)), c('entity.contA', 'entity.contB', 'entity.contC'))


  ## With specified variables
  variables <- new("VariableMetadataList", SimpleList(
    new("VariableMetadata",
      variableSpec = new("VariableSpec", variableId = 'contA', entityId = 'entity'),
      dataType = new("DataType", value = 'NUMBER'),
      dataShape = new("DataShape", value = 'CONTINUOUS')),
    new("VariableMetadata",
      variableClass = new("VariableClass", value = 'native'),
      variableSpec = new("VariableSpec", variableId = 'contB', entityId = 'entity'),
      dataType = new("DataType", value = 'NUMBER'),
      dataShape = new("DataShape", value = 'CONTINUOUS'))
  ))

  result <- correlation(data, 'spearman', variables  = variables, verbose = FALSE)
  # Check data (only sample ids)
  dt <- result@data
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('SampleID'))
  expect_equal(nrow(dt), nSamples)
  # Check stats (all correlation outputs)
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('var1','var2','correlationCoef'))
  expect_equal(nrow(stats), (ncol(testOTU) - 1) * length(variables)) # Should be number of taxa * number of metadata vars
  expect_equal(as.character(unique(stats$var1)), names(testOTU)[2:length(names(testOTU))])
  expect_equal(as.character(unique(stats$var2)), c('entity.contA', 'entity.contB'))

  ## With a date <3
    variables <- new("VariableMetadataList", SimpleList(
    new("VariableMetadata",
      variableSpec = new("VariableSpec", variableId = 'contA', entityId = 'entity'),
      dataType = new("DataType", value = 'NUMBER'),
      dataShape = new("DataShape", value = 'CONTINUOUS')),
    new("VariableMetadata",
      variableClass = new("VariableClass", value = 'native'),
      variableSpec = new("VariableSpec", variableId = 'dateA', entityId = 'entity'),
      dataType = new("DataType", value = 'DATE'),
      dataShape = new("DataShape", value = 'CONTINUOUS'))
  ))

  result <- correlation(data, 'spearman', variables = variables, verbose =FALSE)
  # Check data (only sample ids)
  dt <- result@data
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('SampleID'))
  expect_equal(nrow(dt), nSamples)
  # Check stats (all correlation outputs)
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('var1','var2','correlationCoef'))
  expect_equal(nrow(stats), (ncol(testOTU) - 1) * length(variables)) # Should be number of taxa * number of metadata vars
  expect_equal(as.character(unique(stats$var1)), names(testOTU)[2:length(names(testOTU))])
  expect_equal(as.character(unique(stats$var2)), c('entity.contA', 'entity.dateA'))


})


test_that("correlation returns a ComputeResult with the correct slots" , {

  df <- testOTU
  nSamples <- dim(df)[1]
  sampleMetadata <- data.frame(list(
    "entity.SampleID" = df[["entity.SampleID"]],
    "entity.contA" = rnorm(nSamples),
    "entity.contB" = rnorm(nSamples),
    "entity.contC" = rnorm(nSamples),
    "entity.dateA" = sample(seq(as.Date('1999/01/01'), as.Date('2000/01/01'), by="day"), nSamples)
    ))


  data <- microbiomeComputations::AbundanceData(
              data = df,
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')

  variables <- new("VariableMetadataList", SimpleList(
    new("VariableMetadata",
      variableSpec = new("VariableSpec", variableId = 'contA', entityId = 'entity'),
      dataType = new("DataType", value = 'NUMBER'),
      dataShape = new("DataShape", value = 'CONTINUOUS')),
    new("VariableMetadata",
      variableClass = new("VariableClass", value = 'native'),
      variableSpec = new("VariableSpec", variableId = 'dateA', entityId = 'entity'),
      dataType = new("DataType", value = 'DATE'),
      dataShape = new("DataShape", value = 'CONTINUOUS'))
  ))

  ## Pearson with date and numeric
  result <- correlation(data, 'pearson', variables = variables, verbose = FALSE)
  expect_equal(result@parameters, 'method = pearson')
  expect_equal(result@recordIdColumn, 'entity.SampleID')

  ## With spearman
  result <- correlation(data, 'spearman', variables = variables, verbose = FALSE)
  expect_equal(result@parameters, 'method = spearman')
  expect_equal(result@recordIdColumn, 'entity.SampleID')

})

test_that("correlation fails with improper inputs", {

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

  # Fail when we send in only categorical data
  expect_error(correlation(data, verbose=F))
})
