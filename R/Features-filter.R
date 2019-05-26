##' @title Filter featuress based on their rowData
##'
##' @description
##'
##' The `filterFeatures` methods enables users to filter features
##' based on a variable in their `rowData`. The features matching the
##' filter will be returned as a new object of class `Features`. The
##' filters can be provided as instances of class `AnnotationFilter`
##' (see below) or as formulas.
##'
##' @section Variable filters:
##'
##' The variable filters are filters as defined in the
##' [AnnotationFilter] package. In addition to the pre-defined filter,
##' users can arbitrarily set a field on which to operate. These
##' arbitrary filters operate either on a character variables (as
##' `CharacterVariableFilter` objects) or numerics (as
##' `NumericVariableFilters` objects), which can be created with the
##' `VariableFilter` constructor.
##'
##' @rdname filterFeatures
##'
##' @author Laurent Gatto
##'
##' @name filterFeatures
##'
##' @aliases filterFeatures filterFeatures,Features,formula-method filterFeatures,Features,AnnotationFilter-method CharacterVariableFilter NumericVariableFilter VariableFilter
##'
##' @md
##'
##' @examples
##'
##' ## ----------------------------------------
##' ## Creating character and numberic
##' ## variable filters
##' ## ----------------------------------------
##'
##' VariableFilter(field = "my_var",
##'                value = "value_to_keep",
##'                condition = "==")
##'
##' VariableFilter(field = "my_num_var",
##'                value = 0.05,
##'                condition = "<=")
##'
##' example(Features)
##' 
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
NULL



##' @import AnnotationFilter
##' @exportClass CharacterVariableFilter
##' @rdname filterFeatures
setClass("CharacterVariableFilter", contains = "CharacterFilter")
##' @exportClass NumericVariableFilter
##' @rdname filterFeatures
setClass("NumericVariableFilter", contains = "DoubleFilter")




##' @param field `character(1)` refering to the name of the variable
##'     to apply the filter on.
##' 
##' @param value ‘character()’ or ‘integer()’ value for the
##'     `CharacterVariableFilter` and `NumericVariableFilter` filters
##'     respectively.
##'
##' @param condition ‘character(1)’ defining the condition to be used in
##'     the filter. For ‘NumericVariableFilter’, one of ‘"=="’,
##'     ‘"!="’, ‘">"’, ‘"<"’, ‘">="’ or ‘"<="’. For
##'     ‘CharacterVariableFilter’, one of ‘"=="’, ‘"!="’,
##'     ‘"startsWith"’, ‘"endsWith"’ or ‘"contains"’. Default
##'     condition is ‘"=="’.
##' 
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

##' @param object An instance of class [Features].
##'
##' @param filter Either an instance of class [AnnotationFilter] or a
##'     formula.
##'
##' @param ... Additional parameters. Currently ignored.
##'
##' @exportMethod filterFeatures
##'
##' @rdname filterFeatures
setMethod("filterFeatures",
          c("Features", "AnnotationFilter"),
          function(object, filter, ...) 
              filterFeaturesWithAnnotationFilter(object, filter, ...))

##' @rdname filterFeatures
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


##' @importFrom lazyeval f_eval
filterFeaturesWithFormula <- function(object, filter, ...) {
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      tryCatch(lazyeval::f_eval(filter, data = as.list(x)),
                               error = function(e) rep(FALSE, nrow(x)))
                  })
    object[sel, , ]
}
