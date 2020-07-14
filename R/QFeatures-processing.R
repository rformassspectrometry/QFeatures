##' @title QFeatures processing
##'
##' @description
##'
##' This manual page describes common quantitative proteomics data
##' processing methods using [QFeatures] objects. In the following
##' functions, if `object` is of class `QFeatures`, and optional assay
##' index or name `i` can be specified to define the assay (by name of
##' index) on which to operate.
##'
##' The following functions are currently available:
##'
##' - `logTransform(object, base = 2, i, pc = 0)` log-transforms (with
##'   an optional pseudocount offset) the assay(s).
##'
##' - `normalize(object, method, i)` normalises the assay(s) according
##'   to `method` (see Details).
##'
##' - `scaleTransform(object, center = TRUE, scale = TRUE, i)` applies
##'   [base::scale()] to `SummarizedExperiment` and `QFeatures`
##'   objects.
##'
##' - `sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)`
##'   sweeps out array summaries from `SummarizedExperiment` and
##'   `QFeatures` objects. See [base::sweep()] for details.
##'
##' See the *Processing* vignette for examples.
##'
##' @details
##'
##' The `method` parameter in `normalize` can be one of `"sum"`,
##' `"max"`, `"center.mean"`, `"center.median"`, `"div.mean"`,
##' `"div.median"`, `"diff.meda"`, `"quantiles`", `"quantiles.robust`"
##' or `"vsn"`. The [MsCoreUtils::normalizeMethods()] function returns
##' a vector of available normalisation methods.
##'
##' - For `"sum"` and `"max"`, each feature's intensity is divided by
##'   the maximum or the sum of the feature respectively. These two
##'   methods are applied along the features (rows).
##'
##' - `"center.mean"` and `"center.median"` center the respective
##'   sample (column) intensities by subtracting the respective column
##'   means or medians. `"div.mean"` and `"div.median"` divide by the
##'   column means or medians. These are equivalent to `sweep`ing the
##'   column means (medians) along `MARGIN = 2` with `FUN = "-"` (for
##'   `"center.*"`) or `FUN = "/"` (for `"div.*"`).
##'
##' - `"diff.median"` centers all samples (columns) so that they all
##'   match the grand median by subtracting the respective columns
##'   medians differences to the grand median.
##'
##' - Using `"quantiles"` or `"quantiles.robust"` applies (robust) quantile
##'   normalisation, as implemented in [preprocessCore::normalize.quantiles()]
##'   and [preprocessCore::normalize.quantiles.robust()]. `"vsn"` uses the
##'   [vsn::vsn2()] function.  Note that the latter also glog-transforms the
##'   intensities.  See respective manuals for more details and function
##'   arguments.
##'
##' For further details and examples about normalisation, see
##' [MsCoreUtils::normalize_matrix()].
##'
##' @param object An object of class `QFeatures` or `SummarizedExperiment`.
##'
##' @param x An object of class `QFeatures` or `SummarizedExperiment`
##'     in `sweep`.
##'
##' @param base `numeric(1)` providing the base with respect to which
##'     logarithms are computed. Defaults is 2.
##'
##' @param pc `numeric(1)` with a pseudocount to add to the
##'     quantitative data. Useful when (true) 0 are present in the
##'     data. Default is 0 (no effect).
##'
##' @param center `logical(1)` (default is `TRUE`) value or
##'     numeric-alike vector of length equal to the number of columns
##'     of `object`. See [base::scale()] for details.
##'
##' @param scale `logical(1)` (default is `TRUE`) or a numeric-alike
##'     vector of length equal to the number of columns of
##'     `object`. See [base::scale()] for details.
##'
##' @param method `character(1)` defining the normalisation method to
##'     apply. See Details.
##'
##' @param i A numeric vector or a character vector giving the index or the
##'     name, respectively, of the assay(s) to be processed.
##'
##' @param name A `character(1)` naming the new assay name. Defaults
##'     are `logAssay` for `logTransform`, `scaledAssay` for
##'     `scaleTranform` and `normAssay` for `normalize`.
##'
##' @param MARGIN As in [base::sweep()], a vector of indices giving the
##'     extent(s) of `x` which correspond to `STATS`.
##'
##' @param STATS As in [base::sweep()], the summary statistic which is
##'     to be swept out.
##'
##' @param FUN As in [base::sweep()], the function to be used to carry
##'     out the sweep.
##'
##' @param check.margin As in [base::sweep()], a `logical`.  If `TRUE`
##'     (the default), warn if the length or dimensions of `STATS` do
##'     not match the specified dimensions of `x`.  Set to `FALSE` for
##'     a small speed gain when you know that dimensions match.
##'
##' @param ... Additional parameters passed to inner functions.
##'
##' @aliases logTransform logTransform,SummarizedExperiment-method logTransform,QFeatures-method
##' @aliases scaleTransform scaleTransform,SummarizedExperiment-method scaleTransform,QFeatures-method
##' @aliases normalize normalize,SummarizedExperiment-method normalize,QFeatures-method
##' @aliases sweep sweep,SummarizedExperiment-method sweep,QFeatures-method
##' @aliases normalizeMethods
##'
##' @name QFeatures-processing
##'
##' @rdname QFeatures-processing
##'
##' @examples
##'
##' MsCoreUtils::normalizeMethods()
NULL

## -------------------------------------------------------
##   Transformations
## -------------------------------------------------------


##' @exportMethod logTransform
##' @rdname QFeatures-processing
setMethod("logTransform",
          "SummarizedExperiment",
          function(object, base = 2, pc = 0) {
              assay(object) <- log(assay(object) + pc, base)
              object
          })

##' @rdname QFeatures-processing
setMethod("logTransform",
          "QFeatures",
          function(object, i, name = "logAssay", base = 2, pc = 0) {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")
              if (is.numeric(i)) i <- names(object)[[i]]
              object <- addAssay(object,
                                 logTransform(object[[i]], base, pc),
                                 name)
              addAssayLinkOneToOne(object, from = i, to = name)
          })

##' @exportMethod scaleTransform
##' @rdname QFeatures-processing
setMethod("scaleTransform", "SummarizedExperiment",
          function(object, center = TRUE, scale = TRUE) {
              e <- scale(assay(object), center = center, scale = scale)
              attr(e, "scaled:center") <- NULL
              attr(e, "scaled:scale") <- NULL
              assay(object) <- e
              object
          })

##' @rdname QFeatures-processing
setMethod("scaleTransform", "QFeatures",
          function(object, i, name = "scaledAssay", center = TRUE, scale = TRUE) {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")
              if (is.numeric(i)) i <- names(object)[[i]]
              object <- addAssay(object,
                                 scaleTransform(object[[i]], center, scale),
                                 name)
              addAssayLinkOneToOne(object, from = i, to = name)
          })

## -------------------------------------------------------
##   Normalisation (normalize)
## -------------------------------------------------------

##' @importFrom BiocGenerics normalize
##' @exportMethod normalize
##' @rdname QFeatures-processing
setMethod("normalize", "SummarizedExperiment",
          function(object,
                   method,
                   ...) {
              e <- MsCoreUtils::normalize_matrix(assay(object), method, ...)
              rownames(e) <- rownames(assay(object))
              colnames(e) <- colnames(assay(object))
              assay(object) <- e
              object
          })

## normalise <- normalize

##' @rdname QFeatures-processing
setMethod("normalize", "QFeatures",
          function(object, i, name = "normAssay", method, ...) {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")
              if (is.numeric(i)) i <- names(object)[[i]]
              object <- addAssay(object,
                                 normalize(object[[i]], method, ...),
                                 name)
              addAssayLinkOneToOne(object, from = i, to = name)
          })


## -------------------------------------------------------
##   Sweep
## -------------------------------------------------------

sweepSE <- function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
    e <- base::sweep(assay(x), MARGIN, STATS, FUN, check.margin, ...)
    rownames(e) <- rownames(assay(x))
    colnames(e) <- colnames(assay(x))
    assay(x) <- e
    x
}

##' @exportMethod sweep
##' @rdname QFeatures-processing
setMethod("sweep", "SummarizedExperiment",
          function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
              sweepSE(x, MARGIN, STATS, FUN, check.margin, ...))


##' @rdname QFeatures-processing
setMethod("sweep", "QFeatures",
          function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ..., i, name = "sweptAssay") {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")
              if (is.numeric(i)) i <- names(x)[[i]]
              x <- addAssay(x,
                            sweepSE(x[[i]], MARGIN, STATS, FUN, check.margin, ...),
                            name)
              addAssayLinkOneToOne(x, from = i, to = name)
          })
