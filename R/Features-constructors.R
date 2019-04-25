##' Convert tabular data from a spreadsheet or a `data.frame` into a
##' `Features` object.
##'
##' @title Features from tabular data
##' @param table File or object holding the quantitative data. Can be
##'     either a `character(1)` with the path to a text-based
##'     spreadsheet (comma-separated values by default, but see `...`)
##'     or an object that can be coerced to a `data.frame`.
##' @param ecol A `numeric` indicating the indices of the columns to
##'     be used as expression values. Can also be a `character`
##'     indicating the names of the columns. Caution must be taken if
##'     the column names are composed of special characters like `(`
##'     or `-` that will be converted to a `.` by the `read.csv`
##'     function. If `ecol` does not match, the error message will
##'     dislpay the column names as seen by the `read.csv` function.
##' @param fnames An optional `character(1)` or `numeric(1)`
##'     indicating the column to be used as feature names.
##' @param ... Further arguments that can be passed on to `read.csv`.
##' @return An instance of class `Features`.
##'
##' @author Laurent Gatto
##'
##' @importFrom utils read.csv
##'
##' @seealso The [Features] class for an example on how to manipulate use
##'     `readFeatures` and how to further manipulate the resulting data.
##'
##' @export
readFeatures <- function(table, ecol, fnames, ...)  {
    if (is.data.frame(table)) xx <- table
    else {
        args <- list(...)
        args$file <- table
        if ("rownames" %in% names(args)) {
            if (missing(fnames)) fnames <- args$rownames
            args$rownames <- NULL
        }
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
    attr(assay, "idx") <- 1
    fdata <- DataFrame(xx[, -ecol, drop = FALSE])
    ans <- new("Features",
               assays = SimpleList(assay),
               featureData = SimpleList(fdata),
               colData = DataFrame(row.names = colnames(assay)))
    if (!missing(fnames)) {
        fnames <- fnames[1]
        if (is.numeric(fnames))
            fnames <- colnames(xx)[fnames]
        if (is.na(match(fnames, colnames(xx))))
            stop(fnames, " not found among\n",
                 paste(colnames(xx), paste = ", "))
        featureNames(ans) <- fdata[, fnames]
    }
    if (validObject(ans))
        return(ans)
}


##' @export
##' @rdname Features-class
##' @param assays A `SimpleList` containing the object's assays (typically
##'     matrices).
##' @param featureData A `SimpleList` containing the object's feature metadata
##'     as `DataFrame` instances.
##' @param colData A `DataFrame` with column (sample annotations).
##' @param metadata A `list()` with arbitrary object annotations.
##' @return A new instance of class `Features`.
Features <- function(assays = SimpleList(),
                     featureData = SimpleList(),
                     colData = DataFrame(),
                     metadata = list()) {
    if (length(assays)) {
        nrows <- order(sapply(assays, nrow), decreasing = TRUE)
        idx <- 1
        for (i in nrows) {
            attr(assays[[i]], "idx") <- idx
            idx <- idx + 1
        }
    }
    new("Features",
        assays = assays,
        featureData = featureData,
        colData = colData,
        metadata = metadata)
}
