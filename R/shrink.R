invariant_col<- function(x)
    length(unique(x)) == 1

shrink_col <- function(x) {
    if (invariant_col(x))
        return(unique(x))    
    else 
        res <- SimpleList(x)
    res
}

shrink_cols <- function(x) {
    if (nrow(x) == 1) return(x)
    res <- x[1, ]
    for (i in seq_along(x))
        res[[i]] <- shrink_col(x[[i]])
    res
}


##' A long dataframe can be `shrunk` by mergeing certain rows into a
##' single one. Invariant column, i.e columns that have the the same
##' value along all the rows that need to be merged, are shrunk into a
##' new variables containing that invariant value. Otherwise, the new
##' variables is constructed as a `SimpleList` containing all the
##' original values. The grouping of rows, i.e. the rows that need to
##' be shrunk together as one, is defined by a vector.
##'
##' A shrunken `DataFrame` can be grown back using the `grow`
##' function.
##'
##' It is important to note that if new values are added to the
##' shrunken `DataFrame`, these will be repeated to create invariant
##' values, which might not be relevant.
##'
##' Parallel processing is applied automatically using the
##' `BiocParallel` package. See `bpparam()` for the current
##' parallelisation backend, `registered()` for all registered
##' backends, and use `register(...)` to set a new backend.
##'
##' @title Shrinks and grows a `DataFrame`
##' @param x The `DataFrame` to be shrunk.
##' @param k A ‘vector’ of length `nrow(x)` defining the grouping
##'     based on which the `DataFrame` will be shrunk.
##' @param count `logical(1)` specifying of an additional column
##'     (called by default `.n`) with the tally of rows shrunk into on
##'     new row should be added. Note that if already existing, `.n`
##'     will be silently overwritten.
##' @return
##' @author Laurent Gatto
##' @import S4Vectors
##' @import IRanges
##' @import BiocParallel
##' @importFrom methods as
##' @aliases grow
##' @export
##' @examples
##' library("IRanges")
##' k <- sample(20, 1e5, replace = TRUE)
##' df <- DataFrame(k = k,
##'                 x = round(rnorm(length(k)), 2),
##'                 y = seq_len(length(k)),
##'                 z = sample(LETTERS, length(k), replace = TRUE),
##'                 ir = IRanges(seq_along(k), width = 10),
##'                 r = Rle(sample(5, length(k), replace = TRUE)))
##' df
##'
##' df2 <- shrink(df, df$k, count = TRUE)
##' df2
##'
##' df3 <- grow(df2)
##' df3
##'
##' stopifnot(identical(df[order(df$k), ], df3[, 1:6]))
shrink <- function(x, k, count = FALSE) {
    l <- split(x, k)
    res <- do.call(rbind, bplapply(l, shrink_cols))
    for (i in seq_along(res)) {
        .x <- x[1, i]
        if (is.logical(.x))
            res[[i]] <- as(res[[i]], "LogicalList")
        else if (is.integer(.x))
            res[[i]] <- as(res[[i]], "IntegerList")
        else if (is.double(.x))
            res[[i]] <- as(res[[i]], "NumericList")
        else if (is.character(.x))
            res[[i]] <- as(res[[i]], "CharacterList")
        else if (is.factor(.x))
            res[[i]] <- as(res[[i]], "FactorList")        
        else if (inherits(.x, "Rle"))
            res[[i]] <- as(res[[i]], "RleList")
        else if (inherits(.x, "IRanges"))
            res[[i]] <- as(res[[i]], "IRangesList")
    }
    if (count) 
        res[[".n"]] <- lengths(l)
    res
}


grow_col <- function(x) {
    n <- max(sapply(x, lengths))
    res <- x[rep(1, n), ]
    for (i in seq_along(x)) {
        if (inherits(x[[i]], "List"))
            res[[i]] <- unlist(x[[i]])
    }
    res
}


##' @export
##' @rdname shrink
grow <- function(x) {
    l <- split(x, seq_len(nrow(x)))
    l <- lapply(l, grow_col)
    do.call(rbind, l)
}
