## Desired interface
##
## A function similar to dplyr::filter that will filter rows/features
## based on values in any of the assays rowData. The relevant
## feature(s) are then passed to subsetByFeatue to get the relevant
## features in the other (combined) assays.



##' @import AnnotationFilter
##' @exportClass CharacterVariableFilter
##' @rdname filterFeature
setClass("CharacterVariableFilter", contains = "CharacterFilter")
##' @exportClass NumericVariableFilter
##' @rdname filterFeature
setClass("NumericVariableFilter", contains = "DoubleFilter")


##' @export VariableFilter
##' @rdname filterFeature
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

##' @rdname filterFeature
setGeneric("filterFeature", function(object, filter, ...) standardGeneric("filterFeature"))

setMethod("filterFeature",
          c("Features", "AnnotationFilter")
          function(object, filter, ...) {
              filterFeatureWithAnnotationFilter(object, filter, ...)
          })


setMethod("filterFeature",
          c("Features", "formula")
          function(object, filter, ...) {
              filterFeatureWithFormula(...)
          })


filterFeatureWithAnnotationFilter <- function(object, filter, ...) {
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


filterFeatureWithFormula <- function(object, filter, ...) {
    message("filterFeatureWithFormula")
}

## ## Example
##
## gene <- data.frame(foo = 1:10,
##                    symbol = c(letters[1:9], "b"),
##                    seq_name = paste0("chr", c(1, 4, 4, 8, 1, 2, 5, 3, "X", 4)),
##                    bar = paste0(c("a", "b"), 1:10),
##                    stringsAsFactors = FALSE)
##
## gene
##
## doMatch <- function(x, filter) {
##    do.call(condition(filter), list(x[, field(filter)], value(filter)))
## }
##
## doExtract <- function(x, filter) {
##     x[doMatch(x, filter), ]
## }
##
## f1 <- VariableFilter("foo", 5, ">")
## f2 <- VariableFilter("bar", "a", "startsWith")
## doExtract(gene, f1)
## doExtract(gene, f2)
## ## Also considering
## e <- (~startsWith(markers, "Mito"))
## e <- (~ markers != 'unknown')
## sel <- lapply(experiments(fts2), function(x) lazyeval::f_eval(e, data = as.list(rowData(x))))
## fts2[sel, , ]


