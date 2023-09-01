test_that('Comparator validation works', {

  groupABinList <- veupathUtils::BinList(
                            S4Vectors::SimpleList(
                              c(veupathUtils::Bin(
                                binLabel="binA_a"
                              ))
                            )
                          )
  groupBBinList <- veupathUtils::BinList(
                            S4Vectors::SimpleList(
                              c(veupathUtils::Bin(
                                binLabel="binA_b"
                              ))
                            )
                          )

  
  # Ensure Comparator has a variableId
  expect_error(microbiomeComputations::Comparator(
                          variable = veupathUtils::VariableMetadata(
                            dataShape = veupathUtils::DataShape(value="BINARY")
                          ),
                          groupA = groupABinList,
                          groupB = groupBBinList
  ))

  # Ensure Comparator has no overlap in groupA and groupB
  groupBBinList <- veupathUtils::BinList(
                            S4Vectors::SimpleList(
                              c(veupathUtils::Bin(
                                binLabel="binA_a"
                              ))
                            )
                          )
  expect_error(microbiomeComputations::Comparator(
                          variable = veupathUtils::VariableMetadata(
                            dataShape = veupathUtils::DataShape(value="BINARY")
                          ),
                          groupA = groupABinList,
                          groupB = groupBBinList
  ))

  # Ensure Comparator requires bin starts and ends when variable is continuous
  expect_error(microbiomeComputations::Comparator(
                          variable = veupathUtils::VariableMetadata(
                            dataShape = veupathUtils::DataShape(value="CONTINUOUS")
                          ),
                          groupA = groupABinList,
                          groupB = groupBBinList
  ))

  # Ensure Comparator requires both groupA and groupB
  expect_error(microbiomeComputations::Comparator(
                          variable = veupathUtils::VariableMetadata(
                            dataShape = veupathUtils::DataShape(value="CONTINUOUS")
                          ),
                          groupA = groupABinList
  ))

})


test_that("getGroupLabels returns bin labels", {
  # With a continuous variable
  bin1 <- veupathUtils::Bin(binStart='2', binEnd='3', binLabel="[2, 3)")
  bin2 <- veupathUtils::Bin(binStart='3', binEnd='4', binLabel="[3, 4)")
  bin3 <- veupathUtils::Bin(binStart='4', binEnd='5', binLabel="[4, 5)")
  bin4 <- veupathUtils::Bin(binStart='5', binEnd='6', binLabel="[5, 6)")

  groupABins <- veupathUtils::BinList(S4Vectors::SimpleList(c(bin1, bin2)))
  groupBBins <- veupathUtils::BinList(S4Vectors::SimpleList(c(bin3, bin4)))

  comparatorVariable <- microbiomeComputations::Comparator(
                          variable = veupathUtils::VariableMetadata(
                            variableSpec = VariableSpec(
                              variableId = 'contA',
                              entityId = 'entity'
                            ),
                            dataShape = veupathUtils::DataShape(value="CONTINUOUS")
                          ),
                          groupA = groupABins,
                          groupB = groupBBins
  )

  expect_equal(getGroupLabels(comparatorVariable, "groupA"), c("[2, 3)", "[3, 4)"))
  expect_equal(getGroupLabels(comparatorVariable, "groupB"), c("[4, 5)", "[5, 6)"))
})
