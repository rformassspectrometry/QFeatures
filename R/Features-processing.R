normalise_SE <- function(object, method, ...) {
    if (method == "vsn") {
        e <- Biobase::exprs(vsn::vsn2(assay(object), ...))
    } else if (method == "quantiles") {
        e <- preprocessCore::normalize.quantiles(assay(object), ...)
    } else if (method == "quantiles.robust") {
        e <- preprocessCore::normalize.quantiles.robust(assay(object), ...)
    } else if (method == "center.mean") {
        e <- assay(object)
        center <- colMeans(e, na.rm = TRUE)
        e <- sweep(e, 2L, center, check.margin = FALSE, ...)
    } else if (method == "center.median") {
        e <- assay(object)
        center <- apply(e, 2L, median, na.rm = TRUE)
        e <- sweep(e, 2L, center, check.margin = FALSE, ...)
    } else if (method == "diff.median") {
        e <- assay(object)
        med <- median(as.numeric(e), na.rm = TRUE)
        cmeds <- apply(e, 2L, median, na.rm = TRUE)
        e <- sweep(e, 2L, cmeds - med)
    } else {
        e <- assay(object)
        switch(method,
               max = div <- DelayedArray::rowMaxs(e, na.rm = TRUE),
               sum = div <- rowSums(e, na.rm = TRUE))
        e <- e/div
    }
    rownames(e) <- rownames(assay(object))
    colnames(e) <- colnames(assay(object))
    assay(object) <- e
    object
}


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
##' - `logTransform(object, base = 2, i)`
##'
##' - `normalize(object, method, i)` (or `normalise`) normalises the
##'   assay(s) according to `method` (see Details).
##' 
##' -`scaleTransform(object, center = TRUE, scale = TRUE, i)` applies
##'  [base::scale()] to `SummarizedExperiments` and `Features` objects.
##'
##' @details
##' 
##' The `method` parameter in `normalise` (`normalize`) can be one of
##' `"sum"`, `"max"`, `"quantiles"`, `"center.mean"`,
##' `"center.median"`, `"center.median"`, `"quantiles.robust`" or
##' `"vsn"`.  For `"sum"` and `"max"`, each feature's intensity is
##' divided by the maximum or the sum of the feature
##' respectively. These two methods are applied along the features
##' (rows).
##'
##' `"center.mean"` and `"center.median"` translate the respective
##' sample (column) intensities according to the column mean or
##' median. `"diff.median"` translates all samples (columns) so that
##' they all match the grand median. Using `"quantiles"` or
##' `"quantiles.robust"` applies (robust) quantile normalisation, as
##' implemented in [preprocessCore::normalize.quantiles()] and
##' [preprocessCore::normalize.quantiles.robust()]. `"vsn"` uses the
##' [vsn::vsn2()] function.  Note that the latteralso glog-transforms
##' the intensities.  See respective manuals for more details and
##' function arguments.
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
##' @param i The index or name of the assay to be processed.
##'
##' @aliases logTransform logTransform,SummarizedExperiment-method logTransform,Features-method
##'
##' @aliases scaleTransform scaleTransform,SummarizedExperiment-method scaleTransform,Features-method
##' 
##' @name Features-processing
##'
##' @rdname Features-processing
NULL

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
          function(object, base = 2, i, pc = 0) {
              if (missing(i))
                  i  <-  seq_len(length(object))
              for (ii in i)
                  object[[ii]] <- log(object[[ii]], base, pc)
              object
          })

##' @rdname Features-processing
setMethod("normalize", "SummarizedExperiment",
          function(object,
                   method = c("sum", "max", "center.mean",
                              "center.median", "diff.median",
                              "quantiles", "quantiles.robust", "vsn"),
                   ...)
              normalise_SE(object, match.arg(method), ...))

normalise <- normalize

##' @rdname Features-processing
setMethod("normalize", "Features",
          function(object,
                   method = c("sum", "max", "center.mean",
                              "center.median", "diff.median",
                              "quantiles", "quantiles.robust", "vsn"),
                   i) {
              if (missing(i))
                  i  <-  seq_len(length(object))
              for (ii in i)
                  object[[ii]] <- normalize(object[[ii]], method)

              object
          })


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
          function(object, center = TRUE, scale = TRUE, i) {
              if (missing(i))
                  i  <-  seq_len(length(object))
              for (ii in i)
                  object[[ii]] <- scaleTransform(object[[ii]], center, scale)

              object
          })
