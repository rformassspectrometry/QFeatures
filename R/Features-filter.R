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
##' ## ----------------------------------------
##' ## Filter all features that are associated
##' ## to the Mitochondrion in the markers
##' ## feature variable
##' ## ----------------------------------------
##'
##' ## using the forumla interface, exact mathc
##' filterFeatures(fts2, ~ markers == "Mitochondrion")
##'
##' ## using the forumula intefrace, martial match
##' filterFeatures(fts2, ~startsWith(markers, "Mito"))
##'
##' ## using a user-defined character filter
##' filterFeatures(fts2, VariableFilter("markers", "Mitochondrion"))
##'
##' ## ----------------------------------------
##' ## Filter all features that aren't marked
##' ## as unknown (sub-cellular location) in the
##' ## feature variable
##' ## ----------------------------------------
##'
##' ## using a user-defined character filter
##' filterFeatures(fts2, VariableFilter("markers", "unknown", condition = "!="))
##'
##' ## using the forumula interface
##' filterFeatures(fts2, ~ markers != "unknown")
##'
##' ## ----------------------------------------
##' ## Filter features that have a q-value lower
##' ## or equal to 0.001
##' ## ----------------------------------------
##' 
##' ## using a user-defined numeric filter
##' filterFeatures(fts2, VariableFilter("qValue", 0.001, "<="))
##'
##' ## using the formula interface
##' filterFeatures(fts2, ~ qValue <= 0.001)
##'
##' ## ----------------------------------------
##' ## Negative control - filtering for an
##' ## non-existing markers value or a missing
##' ## feature variable, returning empty results
##' ## ----------------------------------------
##'
##' filterFeatures(fts2, VariableFilter("markers", "not"))
##'
##' filterFeatures(fts2, ~ markers == "not")
##'
##' filterFeatures(fts2, VariableFilter("foo", "bar"))
##'
##' filterFeatures(fts2, ~ foo == "bar")


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
                          do.call(condition(filter),
                                  list(x[, field(filter)],
                                       value(filter)))
                      else
                          rep(FALSE, nrow(x))               
                  })
    object[sel, , ]
}


##' @importFrom lazyeval::f_eval
filterFeaturesWithFormula <- function(object, filter, ...) {
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      tryCatch(lazyeval::f_eval(filter, data = as.list(x)),
                               error = function(e) rep(FALSE, nrow(x)))
                  })
    object[sel, , ]
}
