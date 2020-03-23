##' @title Features processing
##'
##' @description
##'
##' This manual page describes common quantitative proteomics data
##' processing methods using [Features] objects. In the following
##' functions, if `object` is of class `Features`, and optional assay
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
##'   [base::scale()] to `SummarizedExperiments` and `Features` objects.
##'
##' See the *Processing* vignette for examples.
##'
##' @details
##' 
##' The `method` parameter in `normalize` can be one of `"sum"`,
##' `"max"`, `"center.mean"`, `"center.median"`, `"div.mean"`,
##' `"div.median"`, `"diff.meda"`, `"quantiles`", `"quantiles.robust`"
##' or `"vsn"`. The [MsCoreUtils::normalizeMethods()] function
##' returns a vector of available normalisation methods.
##'
##' - For `"sum"` and `"max"`, each feature's intensity is divided by the
##'   maximum or the sum of the feature respectively. These two methods are
##'   applied along the features (rows).
##'
##' - `"center.mean"` and `"center.median"` center the respective sample
##'   (column) intensities by subtracting the respective column means or
##'   medians. `"div.mean"` and `"div.median"` divide by the column means or
##'   medians.
##'
##' - `"diff.median"` centers all samples (columns) so that they all match the
##'   grand median by subtracting the respective columns medians differences to
##'   the grand median.
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
##' @param  object An object of class `Features` or `SummarizedExperiment`.
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
##' @param ... Additional parameters passed to inner functions.
##'
##' @aliases logTransform logTransform,SummarizedExperiment-method logTransform,Features-method
##'
##' @aliases scaleTransform scaleTransform,SummarizedExperiment-method scaleTransform,Features-method
##' 
##' @aliases normalize normalize,SummarizedExperiment-method normalize,Features-method
##'
##' @aliases normalizeMethods
##'
##' @name Features-processing
##'
##' @rdname Features-processing
##'
##' @examples
##'
##' MsCoreUtils::normalizeMethods()
NULL

## -------------------------------------------------------
##   Transformations
## -------------------------------------------------------


##' @exportMethod logTransform
##' @rdname Features-processing
setMethod("logTransform",
          "SummarizedExperiment",
          function(object, base = 2, pc = 0) {
              assay(object) <- log(assay(object) + pc, base)
              object
          })

##' @rdname Features-processing
setMethod("logTransform",
          "Features",
          function(object, i, name = "logAssay", base = 2, pc = 0) {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")  
              addAssay(object,
                       logTransform(object[[i]], base, pc),
                       name)
          })

##' @exportMethod scaleTransform
##' @rdname Features-processing
setMethod("scaleTransform", "SummarizedExperiment",
          function(object, center = TRUE, scale = TRUE) {
              e <- scale(assay(object), center = center, scale = scale)
              attr(e, "scaled:center") <- NULL
              attr(e, "scaled:scale") <- NULL              
              assay(object) <- e
              object
          })

##' @rdname Features-processing
setMethod("scaleTransform", "Features",
          function(object, i, name = "scaledAssay", center = TRUE, scale = TRUE) {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")
              addAssay(object,
                       scaleTransform(object[[i]], center, scale),
                       name)
          })

## -------------------------------------------------------
##   Normalisation (normalize)
## -------------------------------------------------------

##' @importFrom BiocGenerics normalize
##' @exportMethod normalize
##' @rdname Features-processing
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

##' @rdname Features-processing
setMethod("normalize", "Features",
          function(object, i, name = "normAssay", method, ...) {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")
              addAssay(object,
                       normalize(object[[i]], method, ...),
                       name)
          })
