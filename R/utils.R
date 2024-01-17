## Microbiome-specific helpers


## Make default tree when none supplied
## NOT phase 1 - leaving for a later phase
# makeDefaultTree <- function(taxonomy_df) {
  
#     #### Eventually needs to be only the actual variance details
#     #### Eventually needs to make the formula based on col names
#     # phylo.formula <- ...
#     tree <- ape::as.phylo.formula(~Kingdom/Phylum/Class/Order, tax_df)
  
#     nEdges <- dim(tree$edge)[1]
    
#     # Assign edge lengths to all equal 1. 
#     tree$edge.length <- rep(1, nEdges)
    
#     # Root tree using the midpoint
#     # tree <- phytools::pbtree(n=20,scale=30)
#     tree <- phytools::midpoint.root(tree)
    
#     # Force the tree to be ultrametric (all tips equadistant from root)
#     tree <- phytools::force.ultrametric(tree, method="extend")

#     return(tree)
# }

#' @import data.table
#' @import veupathUtils
rankTaxa <- function(df, method=c('median','max','q3','variance')) {

    method <- veupathUtils::matchArg(method)
    #### Notes: Assume df has rows as samples and at least columns Abundance and TaxonomicLevel

    if (identical(method, 'median')) {
      ranked <- df[, list(Abundance=median(Abundance)), by="TaxonomicLevel"]
    } else if (identical(method, 'max')) {
      ranked <- df[, list(Abundance=max(Abundance)), by="TaxonomicLevel"]
    } else if (identical(method, 'q3')) {
      ranked <- df[, list(Abundance=quantile(Abundance, 0.75)), by="TaxonomicLevel"]
    } else if (identical(method, 'variance')) {
      ranked <- df[, list(Abundance=var(Abundance)), by="TaxonomicLevel"]
    } else {
      stop("Unsupported ranking method.")
    }

    data.table::setorderv(ranked, c("Abundance", "TaxonomicLevel"), c(-1, 1))
    
    return(ranked)
}

stripEntityIdFromColumnHeader <- function(columnNames) {
  columnsToFix <- grepl(".", columnNames, fixed=T)
  columnNames[columnsToFix] <- veupathUtils::strSplit(columnNames[columnsToFix], ".", index=2)

  return(columnNames)
}

# Given a data table, a recordIdColumn, and ancestorIdColumns (see slots of AbundanceData or SampleMetadata),
# check to ensure given id columns are valid. Return any errors. 
validateIdColumns <- function(df, record_id_col=character(), ancestor_id_cols=c()) {

  errors <- character()

  if (length(record_id_col) != 1) {
    msg <- "Record ID column must have a single value."
    errors <- c(errors, msg) 
  }

  if (!record_id_col %in% names(df)) {
    msg <- paste("Record ID column is not present in abundance data.frame")
    errors <- c(errors, msg)
  }

  if (!!length(ancestor_id_cols)) {
    if (!all(ancestor_id_cols %in% names(df))) {
      msg <- paste("Not all ancestor ID columns are present in abundance data.frame")
      errors <- c(errors, msg)
    }
  }

  return(errors)
}

isNAorZero <- function(x) {
  return(is.na(x) | x == 0)
}