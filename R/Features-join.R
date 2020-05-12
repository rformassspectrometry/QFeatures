##' This function takes the `SummarizedExperiment`s in `x`, extracts
##' their `rowData`, and returns the column names that are shared and
##' identical between all of them
##'
##' @param x A `Features` of length greater or equal to 2.
##' @return `character()` with the `colData` column names to keep.
##' @author Laurent Gatto
.mcols_to_keep <- function(x) {
    stopifnot(length(x) >= 2)
    mcols <- lapply(experiments(x), rowData)
    mcol_common_cols <- Reduce(intersect, lapply(mcols, names))
    mcols <- lapply(mcols, function(xx) xx[, mcol_common_cols])

    row_names <- lapply(experiments(x), names)
    mcol_common_cols <- Reduce(intersect, row_names)
    
}

joinAssays <- function(x,
                       i,
                       name = "joinedAssay") {
    stopifnot(inherits(x, "Features"))
    if (name %in% names(x))
        stop("Assay of name '", name, "' already exists.")
}
