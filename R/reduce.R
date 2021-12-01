setMethod("showAsCell", "character",
          function (object) {
              n <- 10
              vapply(object, function(x) {
                  if (!is.na(x) & nchar(x) & nchar(x) > n)
                      paste0(substr(x, 1, n), "...")
                  else as.character(x)
              },
              character(1),
              USE.NAMES = FALSE)
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
##' @section Missing values:
##'
##' Missing values do have an important effect on `reduce`. Unless all
##' values to be reduces are missing, they will result in an
##' non-invariant column, and will be dropped with `drop = TRUE`. See
##' the example below.
##'
##' The presence of missing values can have side effects in higher
##' level functions that rely on reduction of `DataFrame` objects.
##'
##' @param x The `DataFrame` to be reduced or expanded.
##'
##' @param k A ‘vector’ of length `nrow(x)` defining the grouping
##'     based on which the `DataFrame` will be shrunk.
##'
##' @param count `logical(1)` specifying of an additional column
##'     (called by default `.n`) with the tally of rows shrunk into on
##'     new row should be added. Note that if already existing, `.n`
##'     will be silently overwritten.
##'
##' @param simplify A `logical(1)` defining if invariant columns
##'     should be converted to simple lists. Default is `TRUE`.
##'
##' @param drop A `logical(1)` specifying whether the non-invariant
##'     columns should be dropped altogether. Default is `FALSE`.
##'
##' @return An expanded (reduced) `DataFrame`.
##'
##' @author Laurent Gatto
##'
##' @import S4Vectors
##'
##' @import IRanges
##'
##' @export reduceDataFrame
##'
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
##'
##' ## Missing values
##' d <- DataFrame(k = rep(1:3, each = 3),
##'                x = letters[1:9],
##'                y = rep(letters[1:3], each = 3),
##'                y2 = rep(letters[1:3], each = 3))
##' d
##'
##' ## y is invariant and can be simplified
##' reduceDataFrame(d, d$k)
##' ## y isn't not dropped
##' reduceDataFrame(d, d$k, drop = TRUE)
##'
##' ## BUT with a missing value
##' d[1, "y"] <- NA
##' d
##'
##' ## y isn't invariant/simplified anymore
##' reduceDataFrame(d, d$k)
##' ## y now gets dropped
##' reduceDataFrame(d, d$k, drop = TRUE)
reduceDataFrame <- function(x, k, count = FALSE,
                            simplify = TRUE, drop = FALSE) {
    res <- split(x, k)
    lens <- unname(lengths(res))
    if (simplify | drop)
        invars <- invariant_cols2(res)
    res <- DataFrame(res)
    if (simplify) {
        ## replace invariant cols
        for (i in invars)
            res[[i]] <- unname(sapply(res[[i]], "[[", 1))
    }
    if (drop)
        res <- res[, invars, drop = FALSE]
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

##' @title Unfold a data frame
##'
##' @description
##'
##' A data frame is said to be *folded* when some cells contain
##' multiple elements. These are often encode as a semi-colon
##' separated character , such as `"a;b"`. This function will
##' transform the data frame to that `"a"` and `"b"` are split and
##' recorded across two lines.
##'
##' The simple example below illustrates a trivial case, where the
##' table below
##'
##' |  X|Y   |
##' |---|----|
##' |  1|a;b |
##' |  2|c   |
##'
##' is unfolded based on the Y variable and becomes
##'
##' |  X|Y  |
##' |---|---|
##' |  1|a  |
##' |  1|b  |
##' |  2|c  |
##'
##' where the value 1 of variable X is now duplicated.
##'
##' If there is a second variable that follows the same pattern as the
##' one used to unfold the table, it also gets unfolded.
##'
##' |  X|Y   | Z  |
##' |---|----|----|
##' |  1|a;b | x;y|
##' |  2|c   | z  |
##'
##' becomes
##'
##' |  X|Y  | Z |
##' |---|---|---|
##' |  1|a  | x |
##' |  1|b  | y |
##' |  2|c  | z |
##'
##' because it is implied that the element in "a;b" are match to "x;y"
##' by their respective indices. Note in the above example, unfolding
##' by Y or Z produces the same result.
##'
##' However, the following table unfolded by Y
##'
##' |  X|Y   |Z   |
##' |--:|:---|:---|
##' |  1|a;b |x;y |
##' |  2|c   |x;y |
##'
##' produces
##'
##' |  X|Y  |Z   |
##' |--:|:--|:---|
##' |  1|a  |x;y |
##' |  1|b  |x;y |
##' |  2|c  |x;y |
##'
##' because "c" and "x;y" along the second row don't match. In this
##' case, unfolding by Z would produce a different result. These
##' examples are also illustrated below.
##'
##' Note that there is no `foldDataFrame()` function. See
##' [reduceDataFrame()] and [expandDataFrame()] to flexibly encode and
##' handle vectors of length > 1 within cells.
##'
##' @param x A `DataFrame` or `data.frame` to be unfolded.
##'
##' @param k `character(1)` referring to a character variable in `x`,
##'     that will be used to unfold `x`.
##'
##' @param split `character(1)` passed to [strsplit()] to split
##'     `x[[k]]`.
##'
##' @return A new object unfolded object of class `class(x)` with
##'     numbers of rows >= `nrow(x)` and columns identical to `x`.
##'
##' @author Laurent Gatto
##'
##' @examples
##'
##' (x0 <- DataFrame(X = 1:2, Y = c("a;b", "c")))
##' unfoldDataFrame(x0, "Y")
##'
##' (x1 <- DataFrame(X = 1:2, Y = c("a;b", "c"), Z = c("x;y", "z")))
##' unfoldDataFrame(x1, "Y")
##' unfoldDataFrame(x1, "Z") ## same
##'
##' (x2 <- DataFrame(X = 1:2, Y = c("a;b", "c"), Z = c("x;y", "x;y")))
##' unfoldDataFrame(x2, "Y")
##' unfoldDataFrame(x2, "Z") ## different
unfoldDataFrame <- function(x, k, split = ";") {
    if (!k %in% names(x))
        stop("Variable '", k, "' not found.")
    k <- k[[1]]
    if (!is.character(x[[k]]))
        stop("'x[[", k, "]]' must be a character.")
    k_name <- k
    ## To unfold by k, first expand/unfold it to a list where each
    ## element of the list contains the unfolded individual elements.
    k <- strsplit(x[[k]], split)

    ## Prepare the result by repeating the rows that need to be
    ## unfolded, i.e. those that were effectively split above.
    kl <- lengths(k)
    i <- rep(seq_len(nrow(x)), kl)
    ans <- x[i, ]

    ## Replace the original variable k by the unfolded version, i.e
    ## c("a;b", "c", "d;e") becomes c("a", "b", "c", "d", "e").
    ans[[k_name]] <- unlist(k)

    ## Replace other character variables if necessary, i.e. they have
    ## the same splitting pattern that k. Otherwise leave them as they
    ## are, i.e. duplicated.
    for (var in setdiff(names(ans), k_name)) {
        if (is.character(ans[[var]])) {
            var_list <- strsplit(x[[var]], split)
            if (identical(lengths(var_list), kl))
                ans[[var]] <- unlist(var_list)
        }
    }
    ans
}
