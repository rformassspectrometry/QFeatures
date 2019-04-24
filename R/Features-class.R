##' @title Quantitative MS Features
##'
##'
##' @aliases Features Features-class class:Features
##'
##' @name Features
##'
##' @description
##'
##' Description comes here ...
##'
##'
##' `Features` objects can be created with the `Features` constructor that takes
##' the following parameters:
##'
##' @details
##'
##' More details come here ...
##'
##' @noRd
##'
##' @import S4Vectors
##' @importFrom SummarizedExperiment assays colData
##' @importFrom Biobase featureNames sampleNames
##'
##' @md
##'
##' @exportClass Features
##'
##' @author Laurent Gatto
##'
##' @examples
##' Features()
setClass("Features",
    representation(
        assays = "SimpleList",
        featreData = "SimpleList",
        colData = "DataFrame",
        metadata = "list"
    ))

setValidity("Features", .valid_Features)

setMethod("show", "Features",
          function(object) {
              selectSome <- S4Vectors:::selectSome
              scat <- function(fmt, vals = character(), exdent = 2, ...) {
                  vals <- ifelse(nzchar(vals), vals, "''")
                  lbls <- paste(S4Vectors:::selectSome(vals), collapse = " ")
                  txt <- sprintf(fmt, length(vals), lbls)
                  cat(strwrap(txt, exdent = exdent, ...), sep="\n")
              }
              if (isEmpty(object)) .show_empty_Features(object)
              else .show_Features(object)
          })


setMethod("names", "Features", function(x) names(x@assays))
setMethod("assays", "Features", function(x, ...) x@assays)
setMethod("dim", "Features", function(x) dim(x@assays[[.main_assay(x)]]))
setMethod("length", "Features", function(x) length(x@assays))
setMethod("isEmpty", "Features", function(x) length(x) == 0)
setMethod("colData", "Features", function(x) x@colData)
setMethod("metadata", "Features", function(x, ...) x@metadata)
setMethod("featureNames", "Features",
          function(object) rownames(object@fData))

featureVars <- function(object, assay = NULL) {
    stopifnot(inherits(object, "Features"))
    if (is.null(assay))
        assay <- .main_assay(object)
    if (is.character(assay)) {
        assay <- assay[1]
        stopifnot(assay %in% names(object))
    }
    names(object@fData[[assay]])
}

setMethod("sampleNames", "Features",
          function(object) rownames(object@colData))
setReplaceMethod("names", "Features",
    function(x, value) {
        names(x@assays) <- value
        x
    })

setReplaceMethod("metadata", "Features",
    function(x, value) {
        if (!is.list(value))
            stop("replacement 'metadata' value must be a list")
        if (!length(value))
            names(value) <- NULL
        x@metadata <- value
        x
    })

## Features <- function(assays = SimpleList(),
##                      fData = DataFrame(),
##                      colData = DataFrame()) {
##     new("Features",
##         assays = assays,
##         fData = fData,
##         colData = colData)
## }
