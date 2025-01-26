##' @rdname QFeatures-class
##'
##' @param colvars A `character()` that selects column(s) in the `colData`.
##'
##' @param rowvars A `character()` with the names of the `rowData`
##'     variables (columns) to retain in any assay.
##'
##' @param index When object is an instance of class `QFeatures`, a `numeric(1)`
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
setMethod("longForm", "QFeatures",
          function(object, colvars = NULL,
                   rowvars = NULL,
                   index = 1L) {
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
                                                i = index)
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
                                            i = index)
              }
          })


## ##' @importFrom reshape2 melt
## setMethod("longForm", "SummarizedExperiment",
##           function(object, colvars = NULL,
##                    rowvars = NULL,
##                    index = NULL) {
##               ## Default is to use all assays in the SE
##               if (is.null(index))
##                   index <- seq_along(assayNames(object))
##               ## Check that index is within bounds
##               if (max(index) > length(assayNames(object)) | min(index) < 1)
##                   stop("Assay index(ces) not within bounds.")
##               if (!is.null(colvars)) {
##                   ## Check that all colvars exist
##                   if (!all(colvars %in% names(colData(object))))
##                       stop("Some 'colvars' not found in colData(.).")
##               }
##               if (!is.null(rowvars)) {
##                   ## Check that all rowvars exist
##                   if (!all(rowvars %in% names(rowData(object))))
##                       stop("Some 'rowvars' not found in rowData(.).")
##               }
##               res <- lapply(index,
##                             function(i) {
##                                 ans <- reshape2::melt(assay(object, i),
##                                                       varnames = c("rowname", "colname"),
##                                                       value.name = "value")
##                                 ans$assayName <- assayNames(object)[i]
##                                 rownames(ans) <- NULL
##                                 ans
##                             })
##               res <- do.call(rbind, res)
##               ## Add colData variables
##               cd <- colData(object)[as.character(res$colname),
##                                     colvars,
##                                     drop = FALSE]
##               rownames(cd) <- NULL
##               res <- cbind(res, cd)
##               ## Add rowData variables
##               rd <- rowData(object)[as.character(res$rowname),
##                                     rowvars,
##                                     drop = FALSE]
##               rownames(rd) <- NULL
##               res <- cbind(res, rd)
##               as(res, "DataFrame")
##           })