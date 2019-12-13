setMethod("showAsCell", "character",
          function (object) {
              n <- 10
              sapply(object, function(x) {
                  if (!is.na(x) & nchar(x) & nchar(x) > n)
                      paste0(substr(x, 1, n), "...")
                  else x
              }, USE.NAMES = FALSE)
          })



## Takes a columns of a CompressedSplitDataFrameList object, iterates
## through its elements and returns FALSE as soon as it finds on
## non-invariant element.
invariant_col2 <- function(x)
    all(lengths(lapply(x, unique)) == 1)    


invariant_cols2 <- function(x) {
    res <- rep(NA, length(x[[1]]))
    for (i in seq_along(res)) {
        res[i] <- invariant_col2(x[, i])
    }
    which(res)
}

##' @title Reduces and expands a `DataFrame`
##'
##' @description
##' 
##' A long dataframe can be *reduced* by mergeing certain rows into a
##' single one.  These new variables are constructed as a `SimpleList`
##' containing all the original values. Invariant columns, i.e columns
##' that have the same value along all the rows that need to be
##' merged, can be shrunk into a new variables containing that
##' invariant value (rather than in list columns). The grouping of
##' rows, i.e. the rows that need to be shrunk together as one, is
##' defined by a vector.
##'
##' The opposite operation is *expand*. But note that for a
##' `DataFrame` to be expanded back, it must not to be simplified.
##'
##' @param x The `DataFrame` to be reduced or expanded.
##' @param k A ‘vector’ of length `nrow(x)` defining the grouping
##'     based on which the `DataFrame` will be shrunk.
##' @param count `logical(1)` specifying of an additional column
##'     (called by default `.n`) with the tally of rows shrunk into on
##'     new row should be added. Note that if already existing, `.n`
##'     will be silently overwritten.
##' @param simplify A `logical(1)` defining if invariant columns
##'     should be converted to simple lists. Default is `TRUE`.
##' @param drop A `logical(1)` specifying whether the non-invariant
##'     columns should be dropped altogether. Default is `FALSE`.
##' @return An expanded (reduced) `DataFrame`.
##' @author Laurent Gatto
##' @import S4Vectors
##' @import IRanges
##' @export reduceDataFrame
##' @examples
##' library("IRanges")
##'
##' k <- sample(100, 1e3, replace = TRUE)
##' df <- DataFrame(k = k,
##'                 x = round(rnorm(length(k)), 2),
##'                 y = seq_len(length(k)),
##'                 z = sample(LETTERS, length(k), replace = TRUE),
##'                 ir = IRanges(seq_along(k), width = 10),
##'                 r = Rle(sample(5, length(k), replace = TRUE)),
##'                 invar = k + 1)
##' df
##'
##' ## Shinks the DataFrame
##' df2 <- reduceDataFrame(df, df$k)
##' df2
##'
##' ## With a tally of the number of members in each group
##' reduceDataFrame(df, df$k, count = TRUE)
##'
##' ## Much faster, but more crowded result
##' df3 <- reduceDataFrame(df, df$k, simplify = FALSE)
##' df3
##'
##' ## Drop all non-invariant columns
##' reduceDataFrame(df, df$k, drop = TRUE)
reduceDataFrame <- function(x, k, count = FALSE, simplify = TRUE, drop = FALSE) {
    res <- split(x, k)
    lens <- lengths(res)
    if (simplify | drop) 
        invars <- invariant_cols2(res)
    res <- DataFrame(res)
    if (simplify) {
        ## replace invariant cols
        for (i in invars)
            res[[i]] <- sapply(res[[i]], "[[", 1)
    }
    if (drop) 
        res <- res[, invars]
    if (count)
        res[[".n"]] <- lens
    res
}


##' @export
##' @rdname reduceDataFrame
expandDataFrame <- function(x, k = NULL) {
    if (is.null(k))
        return(expand(x, recursive = FALSE))
    else
        return(DataFrame(lapply(x, unsplit, k)))
}
