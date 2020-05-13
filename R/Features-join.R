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


##' @export
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
