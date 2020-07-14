##' @exportMethod aggcounts
##' @rdname QFeatures-aggregate
setMethod("aggcounts", "SummarizedExperiment",
          function(object, ...) assay(object, i = "aggcounts", ...))
