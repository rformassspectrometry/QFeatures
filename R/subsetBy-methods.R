## Desired interface
##
## Question: Subset my data `x` by row (protein) `P12345`
##
## Desired output: I want the assays that contain features for that
## protein. These include:
##
## - the `proteins` assay, and specifically the feature names `P12345`.
## - the `peptides` assay, and specifically the peptides associated to
##   protein `P12345`.
## - the `psms` assay, and specifically the PSMs associated to the
##   peptides associated to protein `P12345`.
##
## Desired interface: x["P12345", ]
##
## This requires to
##
## 1. Find the assay that contains feature P12345 (this is the
##    `proteins` assay).
##
## 2. Recursively find the parent assays (these are `peptides` and
##    `psms`).
##
## 3. Find the features in these parent assay that are
##    linked/associated to the protein of interest
##
##
## TODO: Make subsetByFeature work when the feature name matches
##       multiple assays. Probable start by creating a list of feature
##       names, then only subset assays.

find_assay_with_feature_name <- function(x, i) {
    rnms <- rownames(x)
    ans <- lapply(rnms, function(x) i %in% x)
    ans <- sapply(ans, all)
    names(ans)[ans]
}

find_assays_from <- function(x, i) 
    sapply(i, function(ii) names(assayLinks(x, ii)))



setGeneric("subsetByFeature", function(x, y, ...) standardGeneric("subsetByFeature"))

##' @title Subset by feature name
##' 
##' This function will find the assays and features that match
##' directly (by name) or indirectly (through aggregation) the feature
##' name.
##' 
##' The `subsetByFeature` function will first identify the assay that
##' contains the feature(s) `i` and filter the rows matching these
##' feature names exactly. It will then find, in the other assays, the
##' features that produces `i` through aggregation with the
##' `aggregateFeatures` function.
##'
##' See [Features] for an example.
##'
##' @param x An instance of class [Features].
##' @param y A `character` of feature names present in an assay in `x`.
##' @param ... Additional parameters. Ignored.
##' @return An new instance of class [Features] containing relevant
##'     assays and features.
##' @rdname Features-subsetBy
##' 
##' @aliases subsetByFeature subsetByFeature,Features,character-method
##' 
##' @author Laurent Gatto
##' 
##' @exportMethod subsetByFeature
setMethod("subsetByFeature", c("Features", "character"),
          function(x, y, ...) .subsetByFeature(x, y))


.subsetByFeature <- function(x, i) {
    leaf_assay_name  <- find_assay_with_feature_name(x, i)

    if (!length(leaf_assay_name)) 
        stop("Feature not found")

    all_assays_names <- find_assays_from(x, leaf_assay_name)
    all_assays_names <- unique(as.vector(all_assays_names))

    ans <- x[, , all_assays_names]

    ## Let's first collect the feature names for all assays    
    featurename_list <- vector("list", length = length(all_assays_names))
    names(featurename_list) <- all_assays_names

    for (k in leaf_assay_name) 
        featurename_list[[k]] <- i

    for (k in setdiff(all_assays_names, leaf_assay_name)) {
        assay_k <- x[[k]]
        ## which assay(s) created assay_k
        assay_k_parent_name <- names(which(sapply(x@assayLinks, slot, "from") == k))

        for (k2 in assay_k_parent_name) {
            assayLink_k2 <- x@assayLinks[[k2]]@hits
            j <- which(elementMetadata(assayLink_k2)$names_to %in% i)
            featurename_list[[k]] <- union(featurename_list[[k]],
                                           elementMetadata(assayLink_k2)$names_from[j])
        }
        i <- featurename_list[[k]]
    }

    expts <- experiments(x)[featurename_list]
    alnks <- x@assayLinks[all_assays_names]
    alnks <- alnks[featurename_list]

    Features(experiments = expts,
             colData = colData(x),
             sampleMap = sampleMap(x),
             metadata = metadata(x),
             assayLinks = alnks)
}

