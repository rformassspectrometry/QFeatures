##' @description Returns the index of the top n features
##'
##' @param x A `matrix`, typically the assay of a
##'     `SummarizedExperiment` object.
##'
##' @param n A positive `integer(1)` defining the number of top
##'     features to use. Typically 3L.
##'
##' @param fun The `function` that computes a single scalar value for
##'     each row that will be used to order the rows in (by default)
##'     decreasing order (see `decreasing`). Typically `rowSums()`.
##'
##' @param decreasing A `logical(1)` defining whether the values
##'     computed by `fun` need to be ordered in decreasing or
##'     increasing order. Typically `decreasing = TRUE` is used.
##'
##' @param ... Additional parameters passed to `fun`. Typically `na.rm
##'     = TRUE`.
##'
##' @return An `integer` of length `n` times `nlevels(INDEX)`.
##'
##' @noRd
topIdx <- function(x, INDEX, n, fun,
                   decreasing, ...) {
    n <- as.integer(n)[1]
    if (is.na(n) | n < 1)
        stop(sQuote("n"), " must be an integer >= 1.")
    if (nrow(x) != length(INDEX))
        stop(sQuote("nrow(x)"), " and ", sQuote("length(INDEX)"),
             " must be equal.")
    row_summaries <- do.call(fun, list(x, ...))
    o <- order(as.double(row_summaries),
               decreasing = decreasing,
               na.last = TRUE)
    unlist(lapply(split(o, INDEX[o]),
                  head, n),
           use.names = FALSE)
}



##' @param INDEX A `factor` in the sense that `as.factor(INDEX)`
##'     defines how to group the rows of `x`.
##'
##' @param n A positive `integer(1)` defining the number of top
##'     features to use. Typically 3L.
##'
##' @param fun The `function` that computes a single scalar value for
##'     each row that will be used to order the rows in (by default)
##'     decreasing order (see `decreasing`). Default is `rowSums()`.
##'
##' @param decreasing A `logical(1)` defining whether the values
##'     computed by `fun` need to be ordered in decreasing or
##'     increasing order. Default is `decreasing = TRUE`.
##'
##' @param ... Additional parameters passed to `fun`. Typically `na.rm
##'     = TRUE`.
##'
##' @param fcol A `character(1)` naming a rowdata variable (of assay
##'     `i` in case of a `QFeatures`) defining how to group the
##'     features of the assay to select the top `n` ones.
##'
##' @rdname QFeatures-filtering
##'
##' @exportMethod filterTopFeatures
##'
##' @examples
##'
##' ## ----------------------------------------
##' ## Filter top N features: keeps the 2 PSMs
##' ## of each peptides that have the highest
##' ## sum of intensities across all samples
##' ## ----------------------------------------
##'
##' se <- feat1[[1]]
##' filterTopFeatures(se, fcol = "Sequence", n = 2L)
setMethod("filterTopFeatures", "SummarizedExperiment",
          function(object, fcol, n = 3L,
                   fun = rowSums,
                   decreasing = TRUE,
                   ...) {
              if (!fcol %in% names(rowData(object)))
                  stop("'fcol' not found in the assay's rowData.")
              idx <- topIdx(assay(object, 1L),
                            rowData(object)[[fcol]],
                            n, fun, decreasing, ...)
              object[idx, ]
          })

setMethod("filterTopFeatures", "QFeatures",
          function(object, i, fcol, n = 3L,
                   name = "newAssay",
                   fun = rowSums,
                   decreasing = TRUE,
                   ...) {
              if (isEmpty(object))
                  return(object)
              if (name %in% names(object))
                  stop("There's already an assay named '", name, "'.")
              if (missing(i))
                  i <- QFeatures:::main_assay(object)
              ## Create the top-filtered assay
              topAssay <- filterTopFeatures(object[[i]], fcol, n,
                                            fun, decreasing, ...)
              ## Add the top-filtered assay to the QFeatures object
              object <- addAssay(object,
                                 topAssay,
                                 name = name)
              ## Link the input assay to the top-filtered assay
              ## THIS IS NOT CORRECT
              addAssayLink(object,
                           from = i,
                           to  = name,
                           varFrom = fcol,
                           varTo = fcol)
          })
