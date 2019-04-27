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
##' @slot id An `integer(1)` holding the set's identifier. For
##'     internal use at the [Features]-level only.
##'
##' @slot from An `integer(1)` refering to the identifier from which
##'     the current object originates. For internal use at the
##'     [Features]-level only.
##'
##' @slot version A `character(1)` providing the class version. For
##'     internal use only.
##'
##' @seealso [Features] is the main data struture to manipulate and
##'     process quantitative data. 
##' @md
##' 
##' @name FeatureSet
##' @rdname FeatureSet-class
##' @aliases FeatureSet FeatureSet-class class:FeatureSet
##' @exportClass FeatureSet
##'
##' @author Laurent Gatto
##'
##' @examples
##'
##' ## An empty FeatureSet
##' FeatureSet()
##'
##' ## See ?Features for more examples
NULL

setClass("FeatureSet",
         slots = c(assay = "matrix",
                   featureData = "DataFrame",
                   id = "integer",
                   from = "integer",
                   version = "character"),
         prototype = prototype(
             id = NA_integer_,
             from = NA_integer_,
             version = "0.1"))

##' @param assay A `matrix` containing the quantitation data.
##'
##' @param featureData A `DataFrame` with feature annotations.
##'
##' @param id An `integer(1)`. For internal use. 
##' 
##' @export
##' @rdname FeatureSet-class
FeatureSet <- function(assay = matrix(ncol = 0, nrow = 0),
                       featureData = DataFrame()) {
    new("FeatureSet",
        assay = assay,
        featureData = featureData)
}

setMethod("show", "FeatureSet",
          function(object) {
              cat("class:", class(object), "\n")
              cat("dim:", dim(object@assay), "\n")
              scat("Features(%d): %s\n", rownames(object@assay))
              scat("Samples(%d): %s\n", colnames(object@assay))
              scat("Feature variables(%d): %s\n", colnames(object@featureData))
          })

##' @exportMethod sampleNames
##' @rdname FeatureSet-class
##' @param object A `FeatureSet` object.
setMethod("sampleNames", "FeatureSet", function(object) colnames(object@assay))

##' @rdname FeatureSet-class
##' @param value An appropriate replacement value.
setReplaceMethod("sampleNames", c("FeatureSet", "character"),
                 function(object, value) {
                     if (!ncol(object))
                         stop("No samples in empty object")
                     colnames(object@assay) <- value
                     object
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

##' @exportMethod featureVariables
##' @rdname FeatureSet-class
setMethod("featureVariables", c("FeatureSet", "missing"),
          function(x, i, ...) colnames(x@featureData))

##' @exportMethod featureData
##' @rdname FeatureSet-class
setMethod("featureData", c("FeatureSet", "missing"), function(x, i, ...) x@featureData)

##' @exportMethod featureData
##' @rdname FeatureSet-class
setMethod("featureData", c("FeatureSet", "missing"), function(x, i, ...) colnames(x@featureData))

##' @exportMethod featureNames
##' @rdname FeatureSet-class
setMethod("featureNames", c("FeatureSet", "missing"),
          function(x, i, ...) rownames(x@featureData))


##' @exportMethod [
##' @rdname FeatureSet-class
##' @param x The object to subset.
##' @param i Subsetting vector for the object's row.
##' @param j Subsetting vector for the object's columns.
##' @param ... Additional paramaeters (ignored).
##' @param drop Always set to `FALSE`.
setMethod("[", c("FeatureSet", "ANY", "ANY", "missing"),
          function(x, i, j, ..., drop = FALSE) {
              assay2 <- x@assay[i, j, drop = FALSE]
              fd2 <- x@featureData[i, , drop = FALSE]
              FeatureSet(assay = assay2,
                         featureData = fd2,
                         id = x@id)                         
})

.valid_FeatureSet <- function(object) {
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

setValidity("FeatureSet", .valid_FeatureSet)
