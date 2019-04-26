##' Simlple container for an assay and its featre metadata
##'
##' A `FeatureSet` is a simple container to store quantitative data
##' (called the `assay`) and the associated feature metadata (called
##' the `featureData`). The number of rows and the rownames of the
##' assay and feature metadata must match exactly.
##'
##' @slot assay A `matrix` of numerics that stores quantitative data
##'     for a set of features (along the rows) and samples (along the
##'     columns).
##'
##' @slot featureData A `DataFrame` to hold the features annotations.
##'
##' @slot id An `integer(1)` holding the set's identifier or
##'     index. For internal use at the [FeatureList]-level only.
##'
##' @slot version A `character(1)` providing the class version. For
##'     internal use only.
##'
##' @md
##' 
##' @name FeatureSet
##' @rdname FeatureSet-class
##' @aliases FeatureSet FeatureSet-call class:FeatureSet
##'
##' @author Laurent Gatto
##'
##' @examples
##'
##' ## An empty FeatureSet
##' FeatureSet()
##'
##' m <- matrix(1:32, ncol = 4)
##' fd <- DataFrame(A = 1:8, B = letters[1:8])
##' rownames(m) <- rownames(fd) <- paste0("f", 1:8)
##' colnames(m) <- paste0("S", 1:4)
##'
##' ## Example FeatureSets
##' fs1 <- FeatureSet(assay = m,
##'                   featureData = fd)
##' fs1
##'
##' fs2 <- FeatureSet(assay = m[1:2, ],
##'                   featureData = fd[1:2, 1, drop = FALSE])
##' fs2
NULL

setClass("FeatureSet",
         slots = c(assay = "matrix",
                   featureData = "DataFrame",
                   id = "integer",
                   version = "character"),
         prototype = prototype(version = "0.1"))

setMethod("show", "FeatureSet",
          function(object) {
              cat("class:", class(object), "\n")
              cat("dim:", dim(object@assay), "\n")
              scat("Features(%d): %s\n", rownames(object@assay))
              scat("Samples(%d): %s\n", colnames(object@assay))
              scat("Feature variables(%d): %s\n", colnames(object@featureData))
          })


##' @exportMethod dim
##' @rdname FeatureSet-class
setMethod("dim", "FeatureSet", function(x) dim(x@assay))

##' @exportMethod ncol
##' @rdname FeatureSet-class
setMethod("ncol", "FeatureSet", function(x) ncol(x@assay))

##' @exportMethod nrow
##' @rdname FeatureSet-class
setMethod("nrow", "FeatureSet", function(x) nrow(x@assay))

##' @exportMethod assay
##' @rdname FeatureSet-class
setMethod("assay", "FeatureSet", function(x, i, ...) x@assay)

##' @exportMethod featureData
##' @rdname FeatureSet-class
setMethod("featureData", "FeatureSet", function(object) object@featureData)

.valid_FeatureSet_rows <- function(object) {
    rn1 <- rownames(object@assay)
    rn2 <- rownames(object@featureData)
    if (!identical(rn1, rn2))
        stop("Feature names don't match in assay and featureData")

    n1 <- nrow(object@assay)
    n2 <- nrow(object@featureData)
    if (!identical(n1, n2))
        stop("Different number of features in assay and featureData")
    NULL
}

.valid_FeatureSet <- function(object) {
    .valid_FeatureSet_rows(object)
}

setValitidy("FeatureSet", .valid_FeatureSet)


##' @param assay
##'
##' @param featureData
##'
##' @param id An `integer(1)`. For internal use. 
##' 
##' @export
##' @rdname FeatureSet-class
FeatureSet <- function(assay = matrix(ncol = 0, nrow = 0),
                       featureData = DataFrame(),
                       id = NA_integer_) {
    new("FeatureSet",
        assay = assay,
        featureData = featureData,
        id = id[1])
}
                       
