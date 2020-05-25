.merge_2_by_cols <- function(x, y) {
    ## Only keep shared variables
    vars <- intersect(names(x), names(y))
    x <- x[, vars, drop = FALSE]
    y <- y[, vars, drop = FALSE]
    ## Only keep variables that have the same values for matching
    ## columns/rows.
    k <- intersect(rownames(x), rownames(y))
    .x <- x[k, , drop = FALSE]
    .y <- y[k, , drop = FALSE]
    for (j in names(.x)) 
        if (!isTRUE(all.equal(.x[[j]], .y[[j]])))
            x[, j] <- y[, j] <- NULL
    ## Create these to recover rownames
    x$._rownames <- rownames(x)
    y$._rownames <- rownames(y)
    ## Perform full join
    res <- merge(x, y,
                 by = intersect(names(x), names(y)),
                 all.x = TRUE, all.y = TRUE,
                 sort = FALSE)
    ## Set row names and remove temporary column
    rownames(res) <- res[["._rownames"]]
    res[["._rownames"]] <- NULL
    res    
}

##' @importFrom methods as
.merge_2_by_rows <- function(x, y) {
    ## Save class to coerce at the end
    cl <- class(x) 
    res <- merge(x, y,
                 by = 0,
                 all.x = TRUE, all.y = TRUE,
                 sort = FALSE)
    ## Set and remove row names
    rownames(res) <- res[[1]]
    res <- res[, -1]
    as(res, cl)
}

.merge_by_rows <- function(x, y, ...) {
    Reduce(.merge_2_by_rows, list(x, y, ...))
}


.merge_by_cols <- function(x, y, ...) {
    Reduce(.merge_2_by_cols, list(x, y, ...))
}


mergeSElist <- function(x) {
    joined_mcols <- Reduce(.merge_2_by_cols, lapply(x, rowData))
    joined_assay <- Reduce(.merge_2_by_rows, lapply(x, assay))
    joined_coldata <- Reduce(.merge_2_by_cols, lapply(x, colData))
    SummarizedExperiment(joined_assay[rownames(joined_mcols), ],
                         joined_mcols,
                         colData = joined_coldata)
}


##' @title Join assays in a Features object
##'
##' @description
##'
##' This function applies a full-join type of operation on 2 or more
##' assays in a `Features` instance. 
##' 
##' @param x An instance of class [Features].
##' 
##' @param i The indices or names of al least two assays to be joined.
##' 
##' @param name A `character(1)` naming the new assay. Default is
##'     `joinedAssay`. Note that the function will fail if there's
##'     already an assay with `name`.
##' 
##' @return A `Features` object with an additional assay.
##'
##' @details
##'
##' The rows to be joined are chosen based on the rownames of the
##' respective assays. It is the user's responsability to make sure
##' these are meaningful, such as for example refering to unique
##' peptide sequences or proteins. 
##'
##' The join operation acts along the rows and expects the samples
##' (columns) of the assays to be disjoint, i.e. the assays mustn't
##' share any samples. Rows that aren't present in an assay are set to
##' `NA` when merged.
##'
##' The `rowData` slots are also joined. However, only columns that
##' are shared and that have the same values for matching columns/rows
##' are retained. For example of a feature variable `A` in sample `S1`
##' contains value `a1` and variable `A` in sample `S2` in a different
##' assay contains `a2`, then the feature variable `A` is dropped in
##' the merged assay.
##' 
##' @author Laurent Gatto
##' 
##' @export
##'
##' @examples
##' 
##' ## -----------------------------------------------
##' ## An example Features with 3 assays to be joined
##' ## -----------------------------------------------
##' data(feat2)
##' feat2
##'
##' feat2 <- joinAssays(feat2, 1:3)
##'
##' ## Individual assays to be joined, each with 4 samples and a
##' ## variable number of rows.
##' assay(feat2[[1]])
##' assay(feat2[[2]])
##' assay(feat2[[3]])
##' 
##' ## The joined assay contains 14 rows (corresponding to the union
##' ## of those in the initial assays) and 12 samples
##' assay(feat2[["joinedAssay"]])
##'
##' ## The individual rowData to be joined.
##' rowData(feat2[[1]])
##' rowData(feat2[[2]])
##' rowData(feat2[[3]])
##'
##' ## Only the 'Prot' variable is retained because it is shared among
##' ## all assays and the values and coherent across samples (the
##' ## value of 'Prot' for row 'j' is always 'Pj'). The variable 'y' is
##' ## missing in 'assay1' and while variable 'x' is present is all
##' ## assays, the values for the shared rows are different.
##' rowData(feat2[["joinedAssay"]])
joinAssays <- function(x,
                       i,
                       name = "joinedAssay") {
    stopifnot("Object must be of class 'Features'" = inherits(x, "Features"),
              "Need at least 2 assays to join" = length(i) >= 2)
    if (name %in% names(x))
        stop("Assay with name '", name, "' already exists.")
    joined_se <- mergeSElist(as.list(experiments(x[, , i])))
    ## TODO: add the AssayLinks    
    addAssay(x, joined_se, name = name)
}
