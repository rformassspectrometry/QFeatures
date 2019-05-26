##' @title Filter featuress based on their rowData
##'
##' @description
##'
##' The `filterFeatures` methods ...
##'
##' @rdname filterFeatures
##'
##' @author Laurent Gatto
##'
##' @aliases filterFeatures filterFeatures,Features,forumla-method filterFeatures,Features,AnnotationFilter-method
##'
##' @md
##'
##' @examples
##' 
##' filterFeatures(fts2, ~ markers == "Mitochondrion")
##'
##' filterFeatures(fts2, ~startsWith(markers, "Mito"))
##' 
##' filterFeatures(fts2, VariableFilter("markers", "Mitochondrion"))
##'
##' filterFeatures(fts2, VariableFilter("markers", "unknown", condition = "!="))
##'
##' filterFeatures(fts2, ~ markers != "unknown")
##'
##' filterFeatures(fts2, VariableFilter("qValue", 0.001, "<="))
##'
##' filterFeatures(fts2, ~ qValue <= 0.001)


##' @import AnnotationFilter
##' @exportClass CharacterVariableFilter
##' @rdname filterFeatures
setClass("CharacterVariableFilter", contains = "CharacterFilter")
##' @exportClass NumericVariableFilter
##' @rdname filterFeatures
setClass("NumericVariableFilter", contains = "DoubleFilter")


##' @export VariableFilter
##' @rdname filterFeatures
VariableFilter <- function(field,
                           value,
                           condition = "==") {
    if (is.numeric(value))
        new("NumericVariableFilter",
            field = as.character(field),
            value = value,
            condition = condition)                
    else if (is.character(value))
        new("CharacterVariableFilter",
            field = as.character(field),
            value = value,
            condition = condition)
    else
        stop("Value type undefined.")
}

##' @rdname filterFeatures
setGeneric("filterFeatures", function(object, filter, ...) standardGeneric("filterFeatures"))

##' @exportMethod filterFeatures
setMethod("filterFeatures",
          c("Features", "AnnotationFilter"),
          function(object, filter, ...) 
              filterFeaturesWithAnnotationFilter(object, filter, ...))


setMethod("filterFeatures",
          c("Features", "formula"),
          function(object, filter, ...)
              filterFeaturesWithFormula(object, filter, ...))

##' @importFrom BiocGenerics do.call
filterFeaturesWithAnnotationFilter <- function(object, filter, ...) {
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      if (field(filter) %in% names(x))
                          res <- do.call(condition(filter),
                                         list(x[, field(filter)],
                                              value(filter)))
                      else
                          res <- rep(FALSE, nrow(x))               
                  })
    object[sel, , ]
}


##' @importFrom lazyeval::f_eval
filterFeaturesWithFormula <- function(object, filter, ...) {
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      browser()
                      lazyeval::f_eval(filter, data = as.list(x))
                  })
    object[sel, , ]
}
