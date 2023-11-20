# Tests for differential abundance methods

test_that('correlation works with two data tables', {
  nSamples = 200

  testData1 <- data.table(
    "contA1" = rnorm(nSamples),
    "contB1" = rnorm(nSamples),
    "contC1" = rnorm(nSamples)
  )

  testData2 <- data.table(
    "contA2" = rnorm(nSamples),
    "contB2" = rnorm(nSamples),
    "contC2" = rnorm(nSamples)
  )

  corrResult <- correlation(testData1, testData2, 'spearman', verbose=F)
  expect_equal(names(corrResult), c("data1","data2","correlationCoef"))
  expect_equal(nrow(corrResult), ncol(testData1) * ncol(testData2))
  expect_true(!all(is.na(corrResult$correlationCoef)))

  corrResult <- correlation(testData1, testData2, 'pearson', verbose=F)
  expect_equal(names(corrResult), c("data1","data2","correlationCoef"))
  expect_equal(nrow(corrResult), ncol(testData1) * ncol(testData2))
  expect_true(!all(is.na(corrResult$correlationCoef)))

  # Test with NAs
  # Should compute using complete cases. 
  # This to try to catch regression back to the default behavior of cor, which returns results of NA if any value is NA
  testData2$contA2[1] <- NA
  corrResult <- correlation(testData1, testData2, 'pearson', verbose=F)
  expect_equal(names(corrResult), c("data1","data2","correlationCoef"))
  expect_equal(nrow(corrResult), ncol(testData1) * ncol(testData2))
  expect_true(!all(is.na(corrResult$correlationCoef)))

})

test_that('correlation returns an appropriately structured result for abundance data vs metadata', {
  df <- testOTU
  nSamples <- dim(df)[1]
  sampleMetadata <- SampleMetadata(
    data = data.frame(list(
      "entity.SampleID" = df[["entity.SampleID"]],
      "entity.contA" = rnorm(nSamples),
      "entity.contB" = rnorm(nSamples),
      "entity.contC" = rnorm(nSamples)
      # "entity.dateA" = sample(seq(as.Date('1999/01/01'), as.Date('2000/01/01'), by="day"), nSamples)
      )),
    recordIdColumn = "entity.SampleID"
  )


  data <- microbiomeComputations::AbundanceData(
              data = df,
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')
  
  ## All numeric sample variables
  result <- correlation(data, sampleMetadata, 'pearson', verbose = FALSE)
  # Check data (only sample ids)
  dt <- result@data
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('SampleID'))
  expect_equal(nrow(dt), nSamples)
  # Check stats (all correlation outputs)
  statsData <- result@statistics@statistics
  expect_s3_class(statsData, 'data.frame')
  expect_equal(names(statsData), c('data1','data2','correlationCoef'))
  expect_equal(nrow(statsData), (ncol(testOTU) - 1) * length(veupathUtils::findNumericCols(sampleMetadata@data))) # Should be number of taxa * number of metadata vars
  expect_equal(as.character(unique(statsData$data1)), names(testOTU)[2:length(names(testOTU))])
  expect_equal(as.character(unique(statsData$data2)), c('entity.contA', 'entity.contB', 'entity.contC'))
  expect_true(all(!is.na(statsData)))


  ## With method = spearman
  result <- correlation(data, sampleMetadata, 'spearman', verbose = FALSE)
  # Check data (only sample ids)
  dt <- result@data
  expect_s3_class(dt, 'data.table')
  expect_equal(names(dt), c('SampleID'))
  expect_equal(nrow(dt), nSamples)
  # Check stats (all correlation outputs)
  statsData <- result@statistics@statistics
  expect_s3_class(statsData, 'data.frame')
  expect_equal(names(statsData), c('data1','data2','correlationCoef'))
  expect_equal(nrow(statsData), (ncol(testOTU) - 1) * length(veupathUtils::findNumericCols(sampleMetadata@data))) # Should be number of taxa * number of metadata vars
  expect_equal(as.character(unique(statsData$data1)), names(testOTU)[2:length(names(testOTU))])
  expect_equal(as.character(unique(statsData$data2)), c('entity.contA', 'entity.contB', 'entity.contC'))
  expect_true(all(!is.na(statsData)))



  
  # ## With a date <3
  #   variables <- new("VariableMetadataList", SimpleList(
  #   new("VariableMetadata",
  #     variableSpec = new("VariableSpec", variableId = 'contA', entityId = 'entity'),
  #     dataType = new("DataType", value = 'NUMBER'),
  #     dataShape = new("DataShape", value = 'CONTINUOUS')),
  #   new("VariableMetadata",
  #     variableClass = new("VariableClass", value = 'native'),
  #     variableSpec = new("VariableSpec", variableId = 'dateA', entityId = 'entity'),
  #     dataType = new("DataType", value = 'DATE'),
  #     dataShape = new("DataShape", value = 'CONTINUOUS'))
  # ))


})


test_that("correlation returns a ComputeResult with the correct slots" , {

  df <- testOTU
  nSamples <- dim(df)[1]
  sampleMetadata <- SampleMetadata(
    data = data.frame(list(
      "entity.SampleID" = df[["entity.SampleID"]],
      "entity.contA" = rnorm(nSamples),
      "entity.contB" = rnorm(nSamples),
      "entity.contC" = rnorm(nSamples)
      # "entity.dateA" = sample(seq(as.Date('1999/01/01'), as.Date('2000/01/01'), by="day"), nSamples)
      )),
    recordIdColumn = "entity.SampleID"
  )


  data <- microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = 'entity.SampleID')

  ## Pearson
  result <- correlation(data, sampleMetadata, 'pearson', verbose = FALSE)
  expect_equal(result@parameters, 'method = pearson')
  expect_equal(result@recordIdColumn, 'entity.SampleID')
  expect_equal(result@statistics@data1Metadata, 'assay')
  expect_equal(result@statistics@data2Metadata, 'sampleMetadata')

  ## With spearman
  result <- correlation(data, sampleMetadata, 'spearman', verbose = FALSE)
  expect_equal(result@parameters, 'method = spearman')
  expect_equal(result@recordIdColumn, 'entity.SampleID')
  expect_equal(result@statistics@data1Metadata, 'assay')
  expect_equal(result@statistics@data2Metadata, 'sampleMetadata')

})

test_that("correlation fails with improper inputs", {

  df <- testOTU
  counts <- round(df[, -c("entity.SampleID")]*1000) # make into "counts"
  counts[ ,entity.SampleID:= df$entity.SampleID]
  nSamples <- dim(df)[1]
  sampleMetadata <- SampleMetadata(
    data = data.frame(list(
      "entity.SampleID" = df[["entity.SampleID"]],
      "entity.binA" = sample(c("binA_a", "binA_b"), nSamples, replace=T),
      "entity.cat2" = sample(c("cat2_a", "cat2_b"), nSamples, replace=T),
      "entity.cat3" = sample(paste0("cat3_", letters[1:3]), nSamples, replace=T),
      "entity.cat4" = sample(paste0("cat4_", letters[1:4]), nSamples, replace=T)
      )),
    recordIdColumn = "entity.SampleID"
  )


  data <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts,
              # sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')

  data@sampleMetadata <- sampleMetadata

  # Fail when we send in only categorical metadata
  expect_error(correlation(data, sampleMetadata, verbose=F))

  # Fail when sample metadata is missing a sample
  sampleMetadata@data <- sampleMetadata@data[-1, ]
  expect_error(corrleation(data, sampleMetadata, verbose=F))
})

test_that("toJSON works as expected for the CorrelationResult class", {
  
  
  df <- testOTU
  nSamples <- dim(df)[1]
  sampleMetadata <- SampleMetadata(
    data = data.frame(list(
      "entity.SampleID" = df[["entity.SampleID"]],
      "entity.contA" = rnorm(nSamples),
      "entity.contB" = rnorm(nSamples),
      "entity.contC" = rnorm(nSamples)
      )),
    recordIdColumn = "entity.SampleID"
  )


  data <- microbiomeComputations::AbundanceData(
              data = df,
              recordIdColumn = 'entity.SampleID')

  result <- correlation(data, sampleMetadata, 'pearson', verbose = FALSE)
  jsonList <- jsonlite::fromJSON(toJSON(result@statistics))
  expect_equal(names(jsonList), c('data1Metadata', 'data2Metadata', 'statistics'))
  expect_equal(class(jsonList$data1Metadata), "character")
  expect_equal(class(jsonList$data2Metadata), "character")
  expect_equal(names(jsonList$statistics), c('data1', 'data2', 'correlationCoef'))
  expect_equal(unname(unlist(lapply(jsonList$statistics, class))), c('character', 'character', 'character'))
})
