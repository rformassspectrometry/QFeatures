##' Convert tabular data from a spreadsheet or a `data.frame` into a
##' `Features` object.
##'
##' @title Features from tabular data
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
##' @param name An optional `character(1)` to name the [FeatureSet] in
##'     the resulting [Features] object.
##'
##' @return An instance of class [Features].
##'
##' @author Laurent Gatto
##'
##' @importFrom utils read.csv
##'
##' @seealso The [Features] class for an example on how to manipulate use
##'     `readFeatures` and how to further manipulate the resulting data.
##'
##' @importFrom methods new validObject
##'
##' @export
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
    }
    fl <- FeatureSet(assay = assay,
                     featureData = fdata)
    cd <- DataFrame(row.names = colnames(assay))    
    ans <- Features(listData = fl,
                    colData = cd)
    ans <- set_featureSet_ids(ans)
    if (!is.null(name))
        names(ans) <- name[1]
    ans
}


##' @export
##' @importFrom methods extends
##' @rdname Features-class
##' @param ... One or multiple [FeatureSet] instances passed to the
##'     `Features` constructor.
##' @param colData A `DataFrame` with column (sample annotations).
##' @param metadata A `list()` with arbitrary object annotations.
##' @return A new instance of class `Features`.
Features <- function(...,
                     colData = NULL,
                     metadata = list()) {
    if (missing(...))
        return(new("Features"))
    args <- list(...)
    if (length(args) == 1L && extends(class(args[[1L]]), "list")) 
        args <- args[[1L]]
    ## Setting ids
    ids <- seq_len(length(args))
    ids <- ids[order(sapply(args, nrow), decreasing = TRUE)]
    for (i in seq_len(length(args)))
        args[[i]]@id <- ids[i]

    if (is.null(colData))
        colData <- DataFrame(row.names = sampleNames(args[[1]]))
    
    new("Features",
        listData = args,
        colData = colData,
        metadata = metadata)
}
