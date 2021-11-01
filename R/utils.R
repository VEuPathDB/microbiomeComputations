## Microbiome-specific helpers


## Make default tree when none supplied
## NOT phase 1 - leaving for a later phase
# makeDefaultTree <- function(taxonomy_df) {
  
#     #### Eventually needs to be only the actual var details
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
rankTaxa <- function(df, method=c('median','max','q3','var')) {

    method <- veupathUtils::matchArg(method)
    #### Notes: Assume df has rows as samples and at least columns Abundance and TaxonomicLevel

    if (identical(method, 'median')) {
      ranked <- df[, list(Abundance=median(Abundance)), by="TaxonomicLevel"]
    } else if (identical(method, 'max')) {
      ranked <- df[, list(Abundance=max(Abundance)), by="TaxonomicLevel"]
    } else if (identical(method, 'q3')) {
      ranked <- df[, list(Abundance=quantile(Abundance, 0.75)), by="TaxonomicLevel"]
    } else if (identical(method, 'var')) {
      ranked <- df[, list(Abundance=var(Abundance)), by="TaxonomicLevel"]
    } else {
      stop("Unsupported ranking method.")
    }

    data.table::setorderv(ranked, c("Abundance", "TaxonomicLevel"), c(-1, 1))
    
    return(ranked)
}


writeAppResultsToJson <- function(appResults, pattern = NULL, dir = NULL, verbose = c(TRUE, FALSE)) {
  verbose <- matchArg(verbose)
  
  if (is.null(pattern)) pattern <- 'file'
  if (is.null(dir)) dir <- tempdir()
  
  # Simple formatting
  outJson <- getAppJson(appResults)
  
  outFileName <- basename(tempfile(pattern = pattern, tmpdir = dir, fileext = ".json"))
  write(outJson, outFileName)
  veupathUtils::logWithTime(paste('New json file written:', outFileName), verbose)
  
  return(outFileName)
}

#' @importFrom jsonlite unbox
#' @importFrom jsonlite toJSON
getAppJson <- function(appResults) {
  
  outList <- lapply(appResults, function(x) {
    formattedList <- list()
    formattedList$data <- x
    formattedList$metaData <- list()
    formattedList$metaData$computeDetails <- jsonlite::unbox(attr(x, 'computeDetails'))
    formattedList$metaData$defaultRange <- attr(x, 'defaultRange')
    formattedList$metaData$computedAxisLabel <-jsonlite::unbox(attr(x, 'computedAxisLabel'))
    formattedList$metaData$computedVariables <- attr(x, 'computedVariables')
    formattedList$metaData$computedVariableLabels <- attr(x, 'computedVariableLabels')
    formattedList$metaData$pcoaVariance <- attr(x, 'pcoaVariance')
    return(formattedList)
  })
  
  outJson <- jsonlite::toJSON(outList)
  return(outJson)
}

# Proof of principle. Needs additional nice inputs and tmp directories and such.
writeAppDataToFile <- function(appResults, filename="file.tab") {
  
  # make one large dt
  dt <- Reduce(function(...) merge(..., all = TRUE, by='SampleID'), appResults)
  
  # Write to file
  data.table::fwrite(dt, file=filename)
  
  return(filename)
}
