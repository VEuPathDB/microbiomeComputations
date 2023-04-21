# Tests for differential abundance methods
test_that('differentialAbundance returns a correctly formatted data.table', {
  df <- testOTU
  counts <- round(df[, -c("entity.SampleID")]*1000) # make into "counts"
  counts[ ,entity.SampleID:= df$entity.SampleID]
  nSamples <- dim(df)[1]
  sampleMetadata <- data.frame(list(
    "entity.SampleID" = df[["entity.SampleID"]],
    "entity.binA" = sample(c("binA_a", "binA_b"), nSamples, replace=T)
    # "entity.cat2" = sample(c("cat2_a", "cat2_b"), nSamples, replace=T),
    # "entity.cat3" = sample(paste0("cat3_", letters[1:3]), nSamples, replace=T),
    # "entity.cat4" = sample(paste0("cat4_", letters[1:4]), nSamples, replace=T),
    # "entity.contA" = rnorm(nSamples, sd=5)
    ))


  data <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts,
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')
  
  result <- differentialAbundance(data, comparisonVariable = "entity.binA", groupA = NULL, groupB = NULL, method='DESeq', verbose=F)
  dt <- result@data
  expect_equal(names(dt), c('SampleID'))
  expect_s3_class(dt, 'data.table')
  stats <- result@statistics
  expect_s3_class(stats, 'data.frame')
  expect_equal(names(stats), c('log2foldChange','pValue','adjustedPValue','pointID'))
  expect_equal(unname(unlist(lapply(stats, class))), c('numeric','numeric','numeric','character'))

  # result <- differentialAbundance(data, comparisonVariable = "entity.binA", groupA = NULL, groupB = NULL, method='ANCOMBC', verbose=F)
  # dt <- result@data
  # # expect_equal(nrow(dt), nrow(df)) want nrow(dt) = #unique taxa
  # expect_s3_class(dt, 'data.table')
  # expect_equal(names(dt), c('pointID','foldChange','pValue'))
  # expect_equal(unname(unlist(lapply(dt, class))), c('character','numeric','numeric'))
})
# }

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
    "entity.cat4" = sample(paste0("cat4_", letters[1:4]), nSamples, replace=T),
    "entity.contA" = rnorm(nSamples, sd=5)
    ))


  data <- microbiomeComputations::AbsoluteAbundanceData(
              data = counts,
              sampleMetadata = sampleMetadata,
              recordIdColumn = 'entity.SampleID')

  result <- differentialAbundance(data, comparisonVariable = "entity.binA", groupA = NULL, groupB = NULL, method='DESeq', verbose=F)
#   expect_equal(result@parameters, 'method = DESeq')
#   expect_equal(inherits(result@computedVariableMetadata, "VariableMetadataList"), TRUE)
#   expect_equal(result@computedVariableMetadata[[1]]@variableSpec@variableId, 'foldChange')
#   expect_equal(result@computedVariableMetadata[[1]]@variableSpec@entityId, 'entity')
#   expect_equal(result@computedVariableMetadata[[1]]@displayName, 'Fold Change')
})

# test_that("differentialAbundance fails gracefully", {

#   df <- testOTU

#   data <- microbiomeComputations::AbundanceData(
#             data = df,
#             recordIdColumn = 'entity.SampleID')

#   result <- differentialAbundance(data, method='DESeq', verbose=F)
#   dt <- result@data
#   # expect_equal(nrow(dt), nrow(df)) nrow(dt) should equal number of unique taxa
#   expect_s3_class(dt, 'data.table')
#   expect_equal(names(dt), c('pointID','foldChange','pValue'))
#   expect_true(all(is.na(dt[['foldChange']])))
#   expect_equal(typeof(dt$pointID), c('character'))
#   attr <- attributes(dt)
#   expect_equal(typeof(result@parameters), 'character')
#   expect_equal(typeof(result@computedVariableMetadata[[1]]@variableSpec@variableId), 'character')
#   expect_equal(typeof(result@computedVariableMetadata[[1]]@variableSpec@entityId), 'character')
#   expect_equal(typeof(result@computedVariableMetadata[[1]]@dataType@value), 'character')
#   expect_equal(typeof(result@computedVariableMetadata[[1]]@dataShape@value), 'character')
#   expect_equal(typeof(result@computedVariableMetadata[[1]]@displayName), 'character')
#   expect_equal(typeof(result@computedVariableMetadata[[1]]@displayRangeMin), 'double')
#   expect_equal(typeof(result@computedVariableMetadata[[1]]@displayRangeMax), 'double')
  

# }