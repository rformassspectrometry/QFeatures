##' @rdname QFeatures-class
##'
##' @param colvars A `character()` that selects column(s) in the `colData`.
##'
##' @param rowvars A `character()` with the names of the `rowData`
##'     variables (columns) to retain in any assay.
##'
##' @param i When object is an instance of class `QFeatures`, a `numeric(1)`
##'     indicating what assay within each `SummarizedExperiment` object to
##'     return. Default is `1L`. If object is a `SummarizedExperiment`, a
##'     `numeric()` indicating what assays to pull and convert. Default is to
##'     use all assays.
##'
##' @importFrom MultiAssayExperiment longForm
##' @importFrom reshape2 melt
##'
##' @importFrom BiocGenerics longForm
##'
##' @exportMethod longForm
##'
##' @aliases longForm longForm,QFeatures
##' @aliases longFormat
setMethod("longForm", "QFeatures",
          function(object, colvars = NULL,
                   rowvars = NULL,
                   i = 1L)
              longFormQFeatures(object, colvars, rowvars, i))

##' @rdname QFeatures-class
##'
##' @exportMethod longForm
##'
##' @aliases longForm,SummarizedExperiment
setMethod("longForm", "SummarizedExperiment",
          function(object, colvars = NULL,
                   rowvars = NULL,
                   i = seq_along(assays(object)))
              longFormSE(object, colvars, rowvars, i))

longFormQFeatures <- function(object, colvars = NULL,
                              rowvars = NULL,
                              i = 1L) {
    if (length(i) > 1) {
        warning("'i' must be of length 1 - using first element.")
        i <- i[1]
    }
    if (!is.null(rowvars)) {
        rdNames <- rowDataNames(object)
        misNames <- sapply(rdNames,
                           function (x) any(!rowvars %in% x))
        ## Check that all required
        if (any(misNames))
            stop("Some 'rowvars' not found in assay(s): ",
                 paste0(names(misNames)[misNames], collapse = ", "))
        ## Get long format table with quantification values and colvars
        longDataFrame <-
            MultiAssayExperiment::longForm(
                                      as(object, "MultiAssayExperiment"),
                                      colDataCols = colvars,
                                      i = i)
        ## Get the required rowData
        rds <- lapply(rowData(object),
                      function(rd) rd[, rowvars, drop = FALSE])
        rds <- do.call(rbind, rds)
        ## Merge the rowData to the long table
        cbind(longDataFrame,
              rds[as.character(longDataFrame$rowname), , drop = FALSE])
    } else {
        ## If rowvars is null, return the MAE longForm output
        MultiAssayExperiment::longForm(
                                  as(object, "MultiAssayExperiment"),
                                  colDataCols = colvars,
                                  i = i)
    }
}


##' @importFrom reshape2 melt
longFormSE <- function(object, colvars = NULL, rowvars = NULL,
                       i = seq_along(assays(object))) {
    ## Check that indices are within bounds
    if (max(i) > length(assays(object)) | min(i) < 1)
        stop("Argument 'i' out of (assay) bounds.")
    ## Check that all colvars exist
    if (!is.null(colvars)) {
        if (!all(colvars %in% names(colData(object))))
            stop("Some 'colvars' not found in colData(.).")
    }
    ## Check that all rowvars exist
    if (!is.null(rowvars)) {
        if (!all(rowvars %in% names(rowData(object))))
            stop("Some 'rowvars' not found in rowData(.).")
    }
    ## Need names for the assayNames columns. If the object's assays don't have
    ## any names, use the index i set above.
    if (is.null(nms <- assayNames(object)))
        nms <- i
    res <- lapply(seq_along(i),
                  function(ii) {
                      ans <- reshape2::melt(assay(object, i[ii]),
                                            varnames = c("rowname", "colname"),
                                            value.name = "value",
                                            as.is = TRUE)
                      ans$assayName <- nms[ii]
                      rownames(ans) <- NULL
                      ans
                  })
    res <- do.call(rbind, res)
    if (!is.null(colvars)) { ## Add colData variables.
        ## Need object to have colnames
        if (is.null(colnames(object)))
            colnames(object) <- seq_len(ncol(object))
        cd <- colData(object)[as.character(res$colname),
                              colvars,
                              drop = FALSE]
        rownames(cd) <- NULL
        res <- cbind(res, cd)
    }
    if (!is.null(rowvars)) { ## Add rowData variables
        ## Need object to have rownames
        if (is.null(rownames(object)))
            rownames(object) <- seq_len(nrow(object))
        rd <- rowData(object)[as.character(res$rowname),
                              rowvars,
                              drop = FALSE]
        rownames(rd) <- NULL
        res <- cbind(res, rd)
    }
    ## as(res, "DataFrame")
    res
}

##' @noRd
##'
##' @export
longFormat <- function(object,
                       colvars = NULL,
                       rowvars = NULL,
                       i = 1L) {
    .Defunct("longForm")
}