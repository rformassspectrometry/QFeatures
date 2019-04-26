##' Contains a list of `FetureSet` instances
##'
##' This data structure as for sole purpose to store a list of
##' [FeatureSet] instances. All `FeatureSet` instances in a
##' `FeatureList` must however have the same number of columns.
##' Objects of the class can be constructed with the `FeatureList()`
##' constructor. The `[[` and `[` operators work as expected.
##'
##' @seealso [Features] is the main data struture to manipulate and
##'     process quantitative data.
##'
##' @md
##'
##' @name FeatureList
##' @rdname FeatureList-class
##' @aliases FeatureList FeatureList-class class:FeatureList
##' @exportClass FeatureList
##'
##' @author Laurent Gatto
##'
##' @examples
##'
##' ## An empty FeatureList
##' FeatureList()
##'
##' ## See ?Features for more examples
setClass("FeatureList",
         contains = "SimpleList",
         slots = c(version = "character"),
         prototype = prototype(
             elementType = "FeatureSet",
             version = "0.1"))

.valid_FeaturesList_samples <- function(object) {
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
.valid_FeaturesList_indices <- function(object) {
    ids <- sapply(object, slot, "id")
    if (anyNA(ids)) 
        stop("Indices musn't be NA")
    if (anyDuplicated(ids))
        stop("Duplicated indices")
    NULL
}

.valid_FeatureList <- function(object) {
    if (!isEmpty(object)) {
        if (!all(sapply(object, inherits, "FeatureSet")))
            stop("Not all elements are of class 'FeatureSet'")    
        .valid_FeaturesList_samples(object)
        .valid_FeaturesList_indices(object)
    }
    NULL
}

setValidity("FeatureList", .valid_FeatureList)


##' @export
##' @param ... Individual `FeatureSet` instances or a list thereof.
##' @rdname FeatureList-class
##' @importFrom methods extends
FeatureList <- function(...) {
    args <- list(...)
    if (length(args) == 1L && extends(class(args[[1L]]), "list")) 
        args <- args[[1L]]
    ## Setting ids
    ids <- seq_len(length(args))
    ids <- ids[order(sapply(args, nrow), decreasing = TRUE)]
    for (i in seq_len(length(args)))
        args[[i]]@id <- ids[i]
    new("FeatureList", listData = args)
}
