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
    unlist(sapply(i, function(ii) names(assayLinks(x, ii))))


##' @title Subset by feature name
##'
##' @description
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
##'
##' @param y A `character` of feature names present in an assay in `x`.
##'
##' @param ... Additional parameters. Ignored.
##'
##' @return An new instance of class [Features] containing relevant
##'     assays and features.
##'
##' @aliases subsetByFeature,Features,character-method
##'
##' @name subsetByFeature
##'
##' @rdname Features-subsetBy
NULL


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
        ## which assay(s) created assay_k
        assay_k_parent_name <-
            names(which(sapply(x@assayLinks, function(al) any(k %in% al@from))))

        for (k2 in assay_k_parent_name) {
            assayLink_k2 <- x@assayLinks[[k2]]@hits
            if (inherits(assayLink_k2, "List")) 
                assayLink_k2 <- assayLink_k2[[k]]
            l <- featurename_list[[k2]]
            j <- which(elementMetadata(assayLink_k2)$names_to %in% l)
            featurename_list[[k]] <- union(featurename_list[[k]],
                                           elementMetadata(assayLink_k2)$names_from[j])
        }
    }
    
    ## Order the assays in featurename_list to match the assay order in x
    ord <- order(match(names(featurename_list), names(x)))
    featurename_list <- featurename_list[ord]
    ## First subset assays, then subset the features of interest. This is 
    ## suggested by the authors of `MultiAssayExperiment` when x contains 
    ## `SingleCellExperiment` assays. 
    ## Cf https://github.com/waldronlab/MultiAssayExperiment/issues/276
    expts <- subsetByAssay(x, names(featurename_list))
    expts <- experiments(subsetByRow(expts, featurename_list))
    ## First subset the `AssayLink`s from the `AssayLinks`, then subset the 
    ## features of interest.
    alnks <- x@assayLinks[names(featurename_list)]
    alnks <- alnks[featurename_list]

    Features(experiments = expts,
             colData = colData(x),
             sampleMap = sampleMap(x),
             metadata = metadata(x),
             assayLinks = alnks)
}
