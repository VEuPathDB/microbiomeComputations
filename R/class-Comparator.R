
check_comparator <- function(object) {
    errors <- character()

        # # Err if the class of groupA does not match that of groupB
    # if (!identical(class(groupA), class(groupB))) {
    #   stop("groupA and groupB must either both be character vectors or both be BinLists.")
    # }

    # Check that groups exist bc im removing the nice feature of inputting a binary var and letting the method do its thing





    #   if (!isComparisonVarBinned) {

    #     # Do not allow duplicated values
    #     if (!!length(intersect(groupA, groupB))) {
    #       veupathUtils::logWithTime("groupA and groupB cannot share members.", verbose)
    #       stop()
    #     }
    #     # Does each group contain at least one string that matches a value in the comparisonValue column?
    #     if (!any(groupA %in% uniqueComparisonVariableValues)) {
    #       stop("At least one value in groupA must exist as a value in the comparisonValue sampleMetadata column.")
    #     }
    #     if (!any(groupB %in% uniqueComparisonVariableValues)) {
    #       stop("At least one value in groupB must exist as a value in the comparisonValue sampleMetadata column.")
    #     }

    #     # Alert that we discard groupA/B values that are not found in the comparisonVariable
    #     if (!all(groupA %in% uniqueComparisonVariableValues)) {
    #       veupathUtils::logWithTime("Found values in groupA that do not exist in the comparisonVariable. Removing these values.", verbose)
    #       groupA <- groupA[groupA %in% uniqueComparisonVariableValues]
    #     }
    #     if (!all(groupB %in% uniqueComparisonVariableValues)) {
    #       veupathUtils::logWithTime("Found values in groupB that do not exist in the comparisonVariable. Removing these values.", verbose)
    #       groupB <- groupB[groupB %in% uniqueComparisonVariableValues]
    #     }
    #   }

    
    return(if (length(errors) == 0) TRUE else errors)
}

#' Comparator
#' 
#' A class for representing a variable that will be used to compare samples between two groups. The variable's
#' values will be used to split samples into groups.
#' 
#' @slot variable A VariableMetadata
#' @slot groupA BinList
#' @slot groupB BinList
#' @name Comparator-class
#' @rdname Comparator-class
#' @export
Comparator <- setClass("Comparator", representation(
    variable = 'VariableMetadata',
    groupA = 'BinList',
    groupB = 'BinList'
), prototype = prototype(
    variable = new("VariableMetadata"),
    groupA = new("BinList"),
    groupB = new("BinList")
), validity = check_comparator)