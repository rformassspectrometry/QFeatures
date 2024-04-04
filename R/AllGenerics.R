##' @title Placeholder for generics functions documentation
##'
##' @name AllGenerics
##' @rdname AllGenerics
NULL

setGeneric("readQFeatures",
           function(assayData, colAnnotation, ...)
               standardGeneric("readQFeatures"))

setGeneric("subsetByFeature",
           function(x, y, ...) standardGeneric("subsetByFeature"))

setGeneric("zeroIsNA",
           function(object, i) standardGeneric("zeroIsNA"))

setGeneric("infIsNA",
           function(object, i) standardGeneric("infIsNA"))

setGeneric("nNA",
           function(object, i) standardGeneric("nNA"))

setGeneric("logTransform",
           function(object, ...) standardGeneric("logTransform"))

setGeneric("scaleTransform",
           function(object, ...) standardGeneric("scaleTransform"))

setGeneric("aggcounts",
           function(object, ...) standardGeneric("aggcounts"))

## base::sweep
.sweep.useAsDefault <- function(x, MARGIN, STATS, FUN = "-",
                                check.margin = TRUE, ...)
    base::sweep(x, MARGIN, STATS, FUN, check.margin, ...)

setGeneric("sweep",
           signature = "x",
           function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
               standardGeneric("sweep"),
           useAsDefault = .sweep.useAsDefault)
