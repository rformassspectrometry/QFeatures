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


##' This function will find the assays and features that match
##' directly (by name) or indirectly (through aggregation) the feature
##' name.
##' 
##' The `subsetByFeature` function will first identify the assay that
##' contains the feature(s) `i` and filter the rows matching these
##' feature names exactly. It will then find, in the other assays, the
##' features that produces `i` through aggregation with the
##' `combineFeatures` function.
##'
##' See [Features] for an example.
##'
##' @title Subset by feature name
##' @param x An instance of class [Features].
##' @param i Feature names present in one assay in `x`.
##' @return An new instance of class [Features] containing relevant
##'     assays and features.
##' @rdname Features-subsetBy
##' @aliases subsetByFeatures
##' @author Laurent Gatto
##' @export
subsetByFeature <- function(x, i) {
    stopifnot(inherits(x, "Features"))
    stopifnot(is.character(i))
    leaf_assay_name  <- find_assay_with_feature_name(x, i)

    if (!length(leaf_assay_name)) 
        stop("Feature not found")

    if (length(leaf_assay_name) > 1) 
        stop("Feature identified in multiple assays.")

    all_assays_names <- find_assays_from(x, leaf_assay_name)[, 1]

    ans <- x[, , all_assays_names]
    
    ## subset features in leaf_assay (number 1)
    leaf_assay <- x[[leaf_assay_name]]
    leaf_assay <- leaf_assay[i, ]
    ans[[leaf_assay_name]] <- leaf_assay
    leaf_hits <- x@assayLinks[[leaf_assay_name]]@hits
    k <- which(elementMetadata(leaf_hits)$names_to %in% i)
    ans@assayLinks[[leaf_assay_name]]@hits <- leaf_hits[k, ]
    
    ## subset parent assays, iterating over from leaf assay's parents:
    ## number 2:length(ans)
    for (j in 2:length(ans)) {
        assay_j_name <- all_assays_names[j]
        parent_assay_name <- all_assays_names[j-1]
        parent_hits <- x@assayLinks[[parent_assay_name]]@hits
        assay_j <- ans[[j]]

        ## subset assay_j
        k <- which(elementMetadata(parent_hits)$names_to %in% i)        
        assay_j <- assay_j[subjectHits(parent_hits)[k], ]        
        ans[[assay_j_name]] <- assay_j

        ## subset hits_j
        hits_j <- ans@assayLinks[[j]]@hits
        k <- which(elementMetadata(hits_j)$names_to %in% rownames(assay_j))
        ans@assayLinks[[j]]@hits <- hits_j[k, ]

        i <- rownames(assay_j)
    }

    if (validObject(ans))
        return(ans)    
}

