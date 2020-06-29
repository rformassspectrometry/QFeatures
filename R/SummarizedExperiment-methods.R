##' @exportMethod aggcounts
##' @rdname Features-aggregate
setMethod("aggcounts", "SummarizedExperiment",
          function(object, ...) assay(object, i = "aggcounts", ...))

