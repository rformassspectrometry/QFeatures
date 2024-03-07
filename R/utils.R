## tidy_DataFrame_columns <- function(object) {
##     for (nm in names(object)) {
##         x <- object[[nm]]
##         if (inherits(x, "List")) {
##             cat("Variable:", nm, "\n")
##             names(x) <- NULL
##             print(x)
##             cat("Do you want to:\n")
##             cat("  Leave as is      (l)\n")
##             cat("  Drop             (d)\n")
##             cat("  Summarise (mean) (s)\n")
##             k <- scan(n = 1L, what = character())
##             k <- match.arg(k, c("l", "d", "s"))
##             object[[nm]] <- switch(k,
##                                    l = object[[nm]],
##                                    d = NULL,
##                                    s = sapply(x, mean))
##         }
##     }
##     object
## }


## TODO: a programmatic version tidyDataFrame, that takes a vector of
## the form c(fvar = "d", ....) and probably one that drops all List
## columns.
##
## This begs the question whether it is necessary at all to have the
## reduced versions, as one should be able to return to the source
## FeatureSet, that contains all the feature variables.
##
## It might be useful to define some session-wide default rules?


##' Split SummarizedExperiment into an ExperimentList
##'
##' The fonction creates an [ExperimentList] containing
##' [SummarizedExperiment] objects from a [SummarizedExperiment]
##' object (also works with [SingleCellExperiment] objects). `f` is
##' used to split `x`` along the rows (`f`` was a feature variable
##' name) or samples/columns (f was a phenotypic variable name). If f
##' is passed as a factor, its length will be matched to nrow(x) or
##' ncol(x) (in that order) to determine if x will be split along the
##' features (rows) or sample (columns). Hence, the length of f must
##' match exactly to either dimension.
##'
##' This function is not exported and was initially available as
##' scp::.splitSCE().
##'
##' @param x a single [SummarizedExperiment] object
##'
##' @param f a factor or a character of length 1. In the latter case,
##'     `f` will be matched to the row and column data variable names
##'     (in that order). If a match is found, the respective variable
##'     is extracted, converted to a factor if needed.
##' @noRd
.splitSE <- function(x, f) {
    ## Check that f is a factor
    if (is.character(f)) {
        if (length(f) != 1)
            stop("'f' must be of lenght one")
        if (f %in% colnames(rowData(x))) {
            f <- rowData(x)[, f]
        }
        else if (f %in% colnames(colData(x))) {
            f <- colData(x)[, f]
        }
        else {
            stop("'", f, "' not found in rowData or colData")
        }
        if (!is.factor(f))
            f <- factor(f)
    }
    ## Check that the factor matches one of the dimensions
    if (!length(f) %in% dim(x))
        stop("length(f) not compatible with dim(x).")
    if (length(f) == nrow(x)) { ## Split along rows
        xl <- lapply(split(rownames(x), f = f), function(i) x[i, ])
    } else { ## Split along columns
        xl <- lapply(split(colnames(x), f = f), function(i) x[, i])
    }
    ## Convert list to an ExperimentList
    do.call(ExperimentList, xl)
}
