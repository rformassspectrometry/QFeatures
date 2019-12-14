##' @title Placeholder for generics functions documentation
##'
##' @name AllGenerics
##' @rdname AllGenerics
NULL


##' @rdname filterFeatures
setGeneric("filterFeatures", function(object, filter, ...) standardGeneric("filterFeatures"))

##' @rdname subsetByFeature
setGeneric("subsetByFeature", function(x, y, ...) standardGeneric("subsetByFeature"))

##' @rdname missing-data
setGeneric("zeroIsNA", function(object, i) standardGeneric("zeroIsNA"))
