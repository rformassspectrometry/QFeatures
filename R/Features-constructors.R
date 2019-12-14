##' @title Features from tabular data
##'
##' @description
##'
##' Convert tabular data from a spreadsheet or a `data.frame` into a
##' `Features` object.
##'
##' @param table File or object holding the quantitative data. Can be
##'     either a `character(1)` with the path to a text-based
##'     spreadsheet (comma-separated values by default, but see `...`)
##'     or an object that can be coerced to a `data.frame`. It is
##'     advised not to encode characters as factors.
##'
##' @param ecol A `numeric` indicating the indices of the columns to
##'     be used as expression values. Can also be a `character`
##'     indicating the names of the columns. Caution must be taken if
##'     the column names are composed of special characters like `(`
##'     or `-` that will be converted to a `.` by the `read.csv`
##'     function. If `ecol` does not match, the error message will
##'     dislpay the column names as seen by the `read.csv` function.
##'
##' @param fnames An optional `character(1)` or `numeric(1)`
##'     indicating the column to be used as feature names.
##'
##' @param ... Further arguments that can be passed on to `read.csv`
##'     except `stringsAsFactors`, which is always `FALSE`.
##'
##' @param name An `character(1)` to name assay. If not set,
##'     `features` is used.
##'
##' @return An instance of class [Features].
##'
##' @author Laurent Gatto
##'
##' @importFrom utils read.csv
##' @import SummarizedExperiment
##'
##' @seealso The [Features] class for an example on how to manipulate use
##'     `readFeatures` and how to further manipulate the resulting data.
##'
##' @importFrom methods new validObject
##'
##' @md
##' @export
##'
##' @examples
##'
##' ## Load a data.frame with PSM-level data
##' data(hlpsms)
##'
##' ## Create the Features object
##' fts2 <- readFeatures(hlpsms, ecol = 1:10, name = "psms")
##' fts2
readFeatures <- function(table, ecol, fnames, ..., name = NULL)  {
    if (is.data.frame(table)) xx <- table
    else {
        args <- list(...)
        args$file <- table
        if ("rownames" %in% names(args)) {
            if (missing(fnames)) fnames <- args$rownames
            args$rownames <- NULL
        }
        args$stringsAsFactors <- FALSE
        xx <- do.call(read.csv, args)
    }
    if (is.character(ecol)) {
        ecol0 <- ecol
        ecol <- match(ecol0, colnames(xx))
        if (any(is.na(ecol)))
            stop("Column identifiers ",
                 paste(ecol0[is.na(ecol)], collapse = ", "),
                 " not recognised among\n",
                 paste(colnames(xx), paste = ", "))
    }
    assay <- as.matrix(xx[, ecol])
    fdata <- DataFrame(xx[, -ecol, drop = FALSE])

    if (!missing(fnames)) {
        fnames <- fnames[1]
        if (is.numeric(fnames))
            fnames <- colnames(xx)[fnames]
        if (is.na(match(fnames, colnames(xx))))
            stop(fnames, " not found among\n",
                 paste(colnames(xx), paste = ", "))
        rownames(fdata) <- rownames(assay) <- fdata[, fnames]
    } else {
        rownames(fdata) <- rownames(assay) <- seq_len(nrow(assay))
    }
    se <- SummarizedExperiment(assay,
                               rowData = fdata)
    cd <- DataFrame(row.names = colnames(assay))
    if (is.null(name))
        name <- "features"
    el <- structure(list(se), .Names = name[1])
    al <- AssayLinks(AssayLink(name = name[1]))
    ans <- MatchedAssayExperiment(el, colData = cd)
    new("Features",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = al)
}



##' @export
##' @rdname Features-class
##' @param ... See `MultiAssayExperiment` for details.
##' @param assayLinks An optional [AssayLinks] object.
Features <- function(..., assayLinks = NULL) {
    ans <- MatchedAssayExperiment(...)
    if (isEmpty(ans)) assayLinks <- AssayLinks()
    else {
        if (is.null(assayLinks))
            assayLinks <- AssayLinks(names = names(ans))
    }
    new("Features",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = assayLinks)
}



##' @param object An instance of class [Features].
##' @param x A single assay or a *named* list of assays.
##' @param name A `character(1)` naming the single assay (default is
##'     `"newAssay"). Ignored if `x` is a list of assays.
##' @param assayLinks An optional [AssayLinks].
##'
##' @md
##'
##' @rdname Features-class
##'
##' @export
addAssay <- function(object,
                     x,
                     name = "newAssay",
                     assayLinks = AssayLinks(names = name)) {
    stopifnot(inherits(object, "Features"))
    el0 <- object@ExperimentList@listData
    if (is.list(x)) el1 <- x
    else el1 <- structure(list(x), .Names = name[1])
    el <- ExperimentList(c(el0, el1))
    smap <- MultiAssayExperiment:::.sampleMapFromData(colData(object), el)
    if (inherits(assayLinks, "AssayLink"))
        assayLinks <- AssayLinks(assayLinks)
    new("Features",
        ExperimentList = el,
        colData = colData(object),
        sampleMap = smap,
        metadata = metadata(object),
        assayLinks = append(object@assayLinks,
                            assayLinks))
}
