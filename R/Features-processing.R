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
##' @param  object An object of class `Features` or `SummarizedExperiment`.
##'
##' @param base `numeric(1)` providing the base with respect to which
##'     logarithms are computed. Defaults is 2.
##'
##' @param i The index or name of the assay to be processed.
##'
##' @aliases logTransform logTransform,SummarizedExperiment-method logTransform,Features-method
##'
##' @name Features-processing
##'
##' @rdname Features-processing
NULL

##' @exportMethod logTransform
##' @rdname Features-processing
setMethod("logTransform",
          "SummarizedExperiment",
          function(object, base = 2) {
              assay(object) <-  log(assay(object), base)
              object
          })

##' @rdname Features-processing
setMethod("logTransform",
          "Features",
          function(object, base = 2, i) {
              if (missing(i)) {
                  for (i in seq_len(length(object)))
                      object[[i]] <- log(object[[i]], base)
              } else {
                  object[[i]]  <- log(object[[i]], base)
              }
              object
          })
