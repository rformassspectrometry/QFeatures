##' @title Placeholder for generics functions documentation
##'
##' @param x A parameter.
##'
##' @param i Another parameter.
##'
##' @param ... Additional parameters.
##' 
##' 
##' @name AllGenerics
##' @rdname AllGenerics
NULL

##' @exportMethod featureData
##' @rdname AllGenerics
setGeneric("featureData", function(x, i, ...) standardGeneric("featureData"))
##' @exportMethod featureNames
##' @rdname AllGenerics
setGeneric("featureNames", function(x, i, ...) standardGeneric("featureNames"))
##' @exportMethod featureVariables
##' @rdname AllGenerics
setGeneric("featureVariables", function(x, i, ...) standardGeneric("featureVariables"))
##' @exportMethod features
##' @rdname AllGenerics
setGeneric("features", function(x, i, ...) standardGeneric("features"))


##' @exportMethod tidyFeatureData
##' @rdname AllGenerics
setGeneric("tidyFeatureData",
           function(object, i, ...) standardGeneric("tidyFeatureData"))
