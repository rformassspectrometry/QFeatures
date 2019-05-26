##' This function combines the quantitative features of an assay,
##' applying a summarisation function (`fun`) to sets of features as
##' defined by the `fcol` feature variable. The new assay's features
##' will be named based on the unique `fcol` values.
##'
##' @title Combine an assay's quantitative features
##' 
##' @param object An instance of class [Features].
##' 
##' @param i The index or name of the assay which features will be
##'     combined the create the new assay.
##' 
##' @param fcol The feature variable of assay `i` defining how to
##'     summerise the features.
##' 
##' @param name A `character(1)` naming the new assay. Default is
##'     `newAssay`. Note that the function will fail if there's
##'     already an assay with `name`.
##' 
##' @param fun A function used for quantitative feature
##'     aggregation. Default is `median`.
##' 
##' @param ... Additional parameters passed the `fun`.
##' 
##' @return A `Features` object with an additional assay.
##'
##' @md
##'
##' @seealso The `Features` vignette provides an extended example.
##'
##' @export
##'
##' @examples
##'
##' ## Create a Features from a PSM-level data.frame
##' example(readFeatures)
##'
##' ## Combine the PSMs into peptide data
##' fts2 <- combineFeatures(fts2, "psms", "Sequence", name = "peptides")
##' fts2
combineFeatures <- function(object,
                            i,
                            fcol,
                            name = "newAssay",
                            fun = median,
                            ... 
                            ) {
    if (isEmpty(object))
        return(object)
    if (name %in% names(object))
        stop("There's already an assay named '", name, "'.")
    if (missing(i))
        i <- main_assay(object)
    .assay <- assay(object, i)
    .rowdata <- rowData(object[[i]])
    if (missing(fcol))
        stop("fcol require.")
    stopifnot(fcol %in% names(.rowdata))
    groupBy <- .rowdata[[fcol]]

    if (anyNA(.assay)) {
        msg <- paste("Your data contains missing values.",
                     "Please read the relevant section in the",
                     "combineFeatures manual page for details the",
                     "effects of missing values on data aggregation.")
        message(paste(strwrap(msg), collapse = "\n"))
    }

    .combined_assay <- combine_assay(.assay, groupBy, fun, ...)
    .combined_rowdata <- Features::reduceDataFrame(.rowdata, .rowdata[[fcol]],
                                                   simplify = TRUE, drop = TRUE,
                                                   count = TRUE)
    se <- SummarizedExperiment(.combined_assay,
                               rowData = .combined_rowdata[rownames(.combined_assay), ])
    hits <- findMatches(rownames(.combined_assay), groupBy)
    elementMetadata(hits)$names_to <- .rowdata[[fcol]][hits@to]
    elementMetadata(hits)$names_from <- rownames(.assay)[hits@to]

    
    assayLinks <- AssayLink(name = name,
                            from = ifelse(is.character(i), i, names(object)[i]),
                            fcol = fcol,
                            hits = hits)
    object <- addAssay(object, se, name = name, assayLinks = assayLinks)

    if (validObject(object))
        object
}


combine_assay <- function(assay, groupBy, fun, ...) 
    do.call(rbind,
            by(assay, groupBy,
               function(.x) apply(.x, 2, fun, ...)))
