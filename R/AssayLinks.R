##' Links between Assays
##'
##' Links between assays within a [Features] object are handled by the
##' `links` slot. It is composed of a `DataFrame` with three columns,
##' documented below. `AssayLinks` can be created with the
##' `AssayLinks()` and `EmptyAssayLinks()` functions.
##'
##' @param name The name of the assay(s).
##'
##' @param from The name of the parent assay, or `NA_character_`, if
##'     not applicable.
##'
##' @param fcol The feature variable of the parent assay used to
##'     generate the current assay (used in
##'     `combineFeatures`). `NA_character_`, if not applicable.
##'
##' @rdname AssayLinks
##'
##' @name AssayLinks
##' 
##' @md
NULL

##' @rdname AssayLinks
##' @export
EmptyAssayLinks <- function()
    DataFrame(row.names = c("name", "from", "fcol"))


##' @rdname AssayLinks
##' @export
AssayLinks <- function(name, from = NULL, fcol = NULL) {
    if (is.null(from))
        from <- rep(NA_character_, length(name))
    if (is.null(fcol))
        fcol <- rep(NA_character_, length(name))
    DataFrame(name = name,
              from = from,
              fcol = fcol,
              row.names = name)
}

get_assay_link <- function(x, i) {
    j <- x@links$name == i
    if (sum(j) == 1) x@links[i, ]
    else stop("Assay not found.")
}

get_parent_assay_link <- function(x, i) {
    i2 <- get_assay_link(x, i)$from
    get_assay_link(x, i2)
}

addAssayLinks <- function(al1, al2) {
    ans <- rbind(al1, al2)
    if (anyDuplicated(ans$name))
        stop("Found duplicated link names.")
    ans    
}
