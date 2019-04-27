## ------------------
## Private functions
## ------------------

main_assay <- function(object)
    which.max(sapply(object, nrow))

## ----------------------
## id and from management
## ----------------------

##' Returns the [FeatureSet]s `id` slots.
##' 
##' @return The [FeatureSet]s `id` slots.
##' @md
##' @rdname FeatureSet-links
get_featureSet_ids <- function(object) 
    sapply(object, slot, "id")


##' This function sets the [FeatureSet]'s `id` slot to
##' `1:length(object)`. The [FeatureSet] with the larges number of
##' features gets `id` 1, next one 2, ... and the smallest one get
##' `length(object)`.
##' 
##' @param object A Features object.
##' @return A Features object with FeatureSet element ids set to
##'     `1:length(object)`.
##' @md
##' @rdname FeatureSet-links
set_featureSet_ids <- function(object) {
    if (isEmpty(object))
        return(object)
    ln <- length(object)
    new_ids <- seq_len(ln)
    o <- order(dims(object)[1, ], decreasing = TRUE)
    new_ids <- new_ids[o]
    for (i in seq_len(ln))
        object@listData[[i]]@id <- new_ids[i]
    if (validObject(object))
        object
}


## ----------------------------
## Internal validity functions
## ----------------------------

.valid_Features_listData <- function(object) {
    ncols <- sapply(object, ncol)
    if (!all(ncols[1] == ncols))
        stop("Number of samples in feature sets must be equal")
    snms1 <- sampleNames(object[[1]])
    snms <- sapply(object, function(x) all(sampleNames(x) == snms1))
    if (!all(snms))
        stop("FeatureSets have different sample names")
    NULL
}

##' @importFrom methods slot
.valid_Features_indices <- function(object) {
    ids <- sapply(object, slot, "id")
    if (anyNA(ids)) 
        stop("Indices musn't be NA")
    if (anyDuplicated(ids))
        stop("Duplicated indices")
    NULL
}

.valid_Features_elementType <- function(object) {
    if (!isEmpty(object)) {
        if (!all(sapply(object, inherits, "FeatureSet")))
            stop("Not all elements are of class 'FeatureSet'")
    }
    NULL
}

.valid_Features_colData <- function(object) {
    if (isEmpty(object)) {
        if (nrow(object@colData) != 0)
            stop("Samples in colData but none in features")
    } else {
        n1 <- nrow(object@colData)
        n2 <- ncol(object[[1]])
        if (n1 != n2)
            stop("Number of samples in features and colData dont' match")
    }
    NULL
}

.valid_Features <- function(object) {
    .valid_Features_elementType(object)
    .valid_Features_listData(object)
    .valid_Features_indices(object)
    .valid_Features_colData(object)
}

setValidity("Features", .valid_Features)
