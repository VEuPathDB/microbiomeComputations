
check_comparator <- function(object) {

    variable <- object@variable
    groupA <- object@groupA
    groupB <- object@groupB

    errors <- character()

    # Check that the variable has a reasonable variable spec
    if (is.na(variable@variableSpec@variableId)) {
      msg <- "Comparator variable needs a variable id"
      errors <- c(errors, msg)
    }
    
    # Check that groups exist
    if (!length(groupA) || !length(groupA)) {
      msg <- "Both groupA and groupB must be defined"
      errors <- c(errors, msg) 
    }

    if (identical(variable@dataShape@value, "CONTINUOUS")) {
      ## Checks for continuous variables

      # Err if variable is continuous but either group is missing a binStart or binEnd
      if (!all(unlist(lapply(groupA, function(bin) {return(!!length(bin@binStart))})))) {
        msg <- "All groupA bins must have a binStart"
        errors <- c(errors, msg) 
      }
      if (!all(unlist(lapply(groupA, function(bin) {return(!!length(bin@binEnd))})))) {
        msg <- "All groupA bins must have a binEnd"
        errors <- c(errors, msg) 
      }
      if (!all(unlist(lapply(groupB, function(bin) {return(!!length(bin@binStart))})))) {
        msg <- "All groupB bins must have a binStart"
        errors <- c(errors, msg) 
      }
      if (!all(unlist(lapply(groupB, function(bin) {return(!!length(bin@binEnd))})))) {
        msg <- "All groupB bins must have a binEnd"
        errors <- c(errors, msg) 
      }
    } else {
      ## Checks for non-continuous variables
      
      # Ensure no values are duplicated between group A and group B
      groupAValues <- getGroupLabels(object, "groupA")
      groupBValues <- getGroupLabels(object, "groupB")
      
      if (!!length(intersect(groupAValues, groupBValues))) {
        msg <- "groupA and groupB cannot share members"
        errors <- c(errors, msg) 
      }

    }
    
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