##' @title Links between Assays
##'
##' @description
##'
##' Links between assays within a [QFeatures] object are handled by an
##' `AssayLinks` object. It is composed by a list of `AssayLink`
##' instances.
##'
##' @section Constructors:
##'
##' Object can be created with the `AssayLink()` and `AssayLinks()`
##' constructors.
##'
##' @section Methods and functions:
##'
##' - `assayLink(x, i)` accesses the AssayLink at position `i` or with
##'    name `i` in the [QFeatures] object `x`.
##'
##' - `parentAssayLinks(x, i, recursive = FALSE)` accesses the
##'   parent(s) `AssayLinks` or assay with index or name `i`.
##'
##' @section Creating links between assays:
##'
##' - `addAssayLink` takes a parent assay and a child assay contained
##'   in the [QFeatures] object and creates a link given a matching
##'   feature variable in each assay's `rowData`. `addAssayLink` also
##'   allows to link an assay from multiple parent assays (see
##'   Examples).
##' - `addAssayLinkOneToOne` links two assays contained in the
##'   [QFeatures] object. The parent assay and the child assay must
##'   have the same size and contain the same rownames (a different
##'   ordering is allowed). The matching is performed based on the row
##'   names of the assays, instead of a supplied variable name in
##'   `rowData`. Providing multiple parents is not supported.
##'
##' @rdname AssayLinks
##'
##' @name AssayLinks
##'
##' @aliases AssayLinks AssayLink AssayLink-class AssayLinks-class class:AssayLinks class:AssayLink show,AssayLink-method [,AssayLink,character-method [,AssayLinks,character-method
##'
##' @md
##'
##' @examples
##'
##' ##-----------------------------
##' ## Creating an AssayLink object
##' ##-----------------------------
##'
##' al1 <- AssayLink(name = "assay1")
##' al1
##'
##' ##------------------------------
##' ## Creating an AssayLinks object
##' ##------------------------------
##'
##' AssayLinks(al1)
##'
##' al2 <- AssayLinks(names = c("Assay1", "Assay2"))
##' al2
##'
##' ##---------------------------------------
##' ## Adding an AssayLink between two assays
##' ##---------------------------------------
##'
##' ## create a QFeatures object with 2 (identical) assays
##' ## see also '?QFeatures'
##' se <- SummarizedExperiment(matrix(runif(20), ncol = 2,
##'                                   dimnames = list(LETTERS[1:10],
##'                                                   letters[1:2])),
##'                            rowData = DataFrame(ID = 1:10))
##' ft <- QFeatures(list(assay1 = se, assay2 = se))
##'
##' ## assay1 and assay2 are not linked
##' assayLink(ft, "assay2") ## 'from' is NA
##' assayLink(ft, "assay1") ## 'from' is NA
##'
##' ## Suppose assay2 was generated from assay1 and the feature variable
##' ## 'ID' keeps track of the relationship between the two assays
##' ftLinked <- addAssayLink(ft, from = "assay1", to = "assay2",
##'                          varFrom = "ID", varTo = "ID")
##' assayLink(ftLinked, "assay2")
##'
##' ## For one-to-one relationships, you can also use
##' ftLinked <- addAssayLinkOneToOne(ft, from = "assay1", to = "assay2")
##' assayLink(ftLinked, "assay2")
##'
##' ##----------------------------------------
##' ## Adding an AssayLink between more assays
##' ##----------------------------------------
##'
##' ## An assay can also be linked to multiple parent assays
##' ## Create a QFeatures object with 2 parent assays and 1 child assay
##' ft <- QFeatures(list(parent1 = se[1:6, ], parent2 = se[4:10, ], child = se))
##' ft <- addAssayLink(ft, from = c("parent1", "parent2"), to = "child",
##'                    varFrom = c("ID", "ID"), varTo = "ID")
##' assayLink(ft, "child")
##'
NULL

setClassUnion("ListorHits", c("Hits", "List"))

##' @exportClass AssayLink
setClass("AssayLink",
         slots = c(name = "character",
                   from = "character",
                   fcol = "character",
                   hits = "ListorHits"))

##' @exportClass AssayLinks
setClass("AssayLinks",
         contains = "SimpleList",
         prototype = prototype(
             elementType = "AssayLink"))

##' @rdname AssayLinks
##' @param object An `AssayLink` object to show.
##' @export
setMethod("show", "AssayLink",
          function(object) {
              cat("AssayLink for assay <", object@name, ">\n",
                  "[from:", paste(object@from, collapse = ","),
                  "|fcol:", paste(object@fcol, collapse = ","),
                  "|hits:", ifelse(inherits(object@hits, "List"),
                                   paste(vapply(object@hits,
                                                length,
                                                numeric(1)),
                                         collapse = ","),
                                   length(object@hits)),
                  "]\n", sep = "")
          })

## --------------
## Constructors
## --------------

##' @param name A mandatory name of the assay(s).
##'
##' @param from The name of the parent assay, or `NA_character_`, if
##'     not applicable.
##'
##' @param fcol The feature variable of the parent assay used to
##'     generate the current assay (used in
##'     `aggregateFeatures`). `NA_character_`, if not applicable.
##'
##' @param hits An object of class [S4Vectors::Hits] matching the
##'     features of two assays.
##'
##' @rdname AssayLinks
##'
##' @md
##'
##' @export
AssayLink <- function(name, from = NA_character_,
                      fcol = NA_character_,
                      hits = Hits()) {
    new("AssayLink", name = name[1],
        from = from, fcol = fcol,
        hits = hits)
}


##' @param ... A set of `AssayLink` objects or a list thereof.
##'
##' @param names A `character()` of `AssayLink` names. If provided,
##'     `...` are ignored, and `names` is used to create an
##'     `AssayLinks` object with `AssayLink` instances with names
##'     `names`.
##'
##' @rdname AssayLinks
##'
##' @importFrom methods extends
##'
##' @md
##'
##' @export
AssayLinks <- function(..., names = NULL) {
    if (!is.null(names))
        return(AssayLinks(sapply(names, AssayLink)))
    args <- list(...)
    if (length(args) == 1L && extends(class(args[[1L]]), "list"))
        args <- args[[1L]]
    names(args) <- vapply(args, slot, "name", FUN.VALUE = character(1))
    new("AssayLinks", listData = args)
}


## -------------------------
## Functions for navigation
## -------------------------


##' @export
##' @rdname AssayLinks
##' @return `assayLink` returns an instance of class `AssayLink`.
assayLink <- function(x, i)
    x@assayLinks[[i]]


##' @param x An instance of class [QFeatures].
##'
##' @param i The index or name of the assay whose `AssayLink` and
##'     parents `AssayLink` instances are to be returned. For `[`, the
##'     feature names to filter on.
##'
##' @return `assayLinks` returns an instance of class `AssayLinks`.
##'
##' @export
##' @md
##' @rdname AssayLinks
assayLinks <- function(x, i) {
    this <- assayLink(x, i)
    if (all(is.na(this@from)))
        return(AssayLinks(this))
    ans <- list()
    froms <- this@from
    while (length(froms)) {
        ans <- append(ans, this)
        this <- assayLink(x, froms[1])
        froms <- c(froms[-1], this@from[!is.na(this@from)])
    }
    ans <- append(ans, this)
    return(AssayLinks(ans))
}

##' @importFrom methods is
##' @param j ignored.
##' @param drop ignored.
##' @rdname AssayLinks
setMethod("[", c("AssayLink", "character"),
          function(x, i, j, ..., drop = TRUE) {
              if (is(x@hits, "List")) {
                  hits <- lapply(x@hits, function(hit) {
                      k <- which(elementMetadata(hit)$names_to %in% i)
                      hit[k, ]
                  })
                  x@hits <- List(hits)
              } else {
                  k <- which(elementMetadata(x@hits)$names_to %in% i)
                  x@hits <-  x@hits[k, ]
              }
              x
          })

##' @rdname AssayLinks
setMethod("[", c("AssayLinks", "list"),
          function(x, i, j, ..., drop = TRUE) {
              stopifnot(identical(names(x), names(i)))
              for (j in names(x)) {
                  alnk <- x[[j]]
                  fnms <- i[[j]]
                  x[[j]] <- alnk[fnms]
              }
              x
          })


## -----------------------------------
## Functions for creating custom links
## -----------------------------------

## The function takes the rowData of an assays to link from, the
## rowData of an assay to link to, the corresponding feature variable
## names in both rowData that link 2 assays together. The function
## returns a `Hits` object.  the corresponding feature variables.
.get_Hits <- function(rdFrom,
                      rdTo,
                      varFrom,
                      varTo) {
    ## Get the shared feature variable
    matchFrom <- unlist(rdFrom[varFrom], use.names = FALSE)
    matchTo <- unlist(rdTo[varTo], use.names = FALSE)
    ## Find hits
    hits <- findMatches(matchFrom, matchTo)
    if (length(c(hits@from, hits@to)) == 0) stop(paste0("No match found"))
    ## Add the row names corresponding to the hits
    elementMetadata(hits)$names_from <- rownames(rdFrom)[hits@from]
    elementMetadata(hits)$names_to <- rownames(rdTo)[hits@to]
    ## Return the Hits object
    return(hits)
}


## The function takes a QFeatures object, the assay name*s* to link from, the
## assay name to link to and the feature variable names contained in the rowData
## of the assay*s* to link from and the assay to link to. The function returns
## an AssayLink object that links the parent assay*s* to the child assay given
## the relationship between the corresponding feature variables.
.create_assay_link <- function(object,
                               from,
                               to,
                               varFrom,
                               varTo) {
    if (any(to %in% from))
        stop("Adding an AssayLink between an assay and itself is not allowed.")
    if (missing(varFrom) | missing(varTo)) {
        ## Create the list of hits between the child assay and the
        ## parent assays based on the rownames
        rdTo <- cbind(rowData(object[[to]]), ._rownames = rownames(object[[to]]))
        hits <- lapply(seq_along(from), function(ii) {
            rdFrom <- cbind(rowData(object[[from[ii]]]),
                            ._rownames = rownames(object[[from[[ii]]]]))
            .get_Hits(rdFrom = rdFrom,
                      rdTo = rdTo,
                      "._rownames", "._rownames")
        })
        varFrom <- "._rownames"
    } else if (length(from) == length(varFrom)) {
        ## Create the list of hits between the child assay and the
        ## parent assays based on the supplied varFrom and varTo
        hits <- lapply(seq_along(from), function(ii) {
            .get_Hits(rdFrom = rowData(object[[from[ii]]]),
                      rdTo = rowData(object[[to]]),
                      varFrom[[ii]], varTo)
        })
    } else
        stop("'from' and 'varFrom' must have same length.")

    ## Format the hits slot to the expected class ("ListorHits")
    if (length(hits) > 1) {
        hits <- List(hits)
        names(hits) <- from
    } else {
        hits <- hits[[1]]
    }

    ## Return a multi-parent AssayLink
    AssayLink(name = to,
              from = from,
              fcol = varFrom,
              hits = hits)
}

## Function that updates the QFeatures object's AssayLinks with the
## provided AssayLink object.
.update_assay_links <- function (object, al) {
    if (!inherits(al ,"AssayLink"))
        stop("'al' must be an AssayLink object.")
    ## Get the hits slot
    hits <- al@hits
    if (inherits(hits, "Hits")) hits <- List(hits)
    ## Check the child indexing on rownames
    isCorrectToLink <- vapply(hits,
                              function(l) all(elementMetadata(l)$names_to %in% rownames(object[[al@name]])),
                              logical(1))
    if (any(!isCorrectToLink))
        stop("Invalid AssayLink. At least one of the 'hits' metadata 'names_to' does not match the rownames.")
    ## Check the parent indexing on rownames
    isCorrectFromLink <- vapply(seq_along(hits),
                                function(i) all(elementMetadata(hits[[i]])$names_from %in% rownames(object[[al@from[i]]])),
                                logical(1))
    if (any(!isCorrectFromLink))
        stop("Invalid AssayLink. The AssayLink metadata 'names_from' does not match the rownames.")

    ## Overwrite the AssayLink
    object@assayLinks@listData[[al@name]] <- al

    if (validObject(object))
        return(object)
}

##' @rdname AssayLinks
##'
##' @param from A `character()` or `integer()` indicating which
##'     assay(s) to link from in `object`
##' @param to A `character(1)` or `integer(1)` indicating which assay
##'     to link to in `object`
##' @param varFrom A `character()` indicating the feature variable(s)
##'     to use to match the `from` assay(s) to the `to`
##'     assay. `varFrom` must have the same length as `from` and is
##'     assumed to be ordered as `from`.
##' @param varTo A `character(1)` indicating the feature variable to
##'     use to match the `to` assay to the `from` assay(s).
##'
##' @export
addAssayLink <- function(object,
                         from,
                         to,
                         varFrom,
                         varTo) {
    if (is.numeric(from)) from <- names(object)[from]
    if (is.numeric(to)) to <- names(object)[[to]]
    ## Create the assay link
    al <- .create_assay_link(object, from, to, varFrom, varTo)
    ## Update the assay link in the QFeatures object
    .update_assay_links(object, al)
}

##' @rdname AssayLinks
##'
##' @export
addAssayLinkOneToOne <- function(object,
                                 from,
                                 to) {
    if (is.numeric(from)) from <- names(object)[[from]]
    if (is.numeric(to)) to <- names(object)[[to]]
    if (length(from) > 1)
        stop("One to one links are not supported for multiple parents.")
    if (any(to %in% from))
        stop("Adding an AssayLink between an assay and itself is not allowed.")

    ## Check that assays have same size
    N <- unique(dims(object)[1, c(from, to)])
    if (length(N) != 1)
        stop("The 'from' and 'to' assays must have the same number of rows.")

    ## Check both assays contain the same rownames (different order is allowed)
    rdFrom <- rowData(object[[from]])
    rdTo <- rowData(object[[to]])
    if (length(intersect(rownames(rdFrom), rownames(rdTo))) != N)
        stop(paste0("Different rownames found in assay '", from,
                    "' and assay '", to, "'."))

    ## Create the linking variable
    rdFrom$._oneToOne <- rownames(rdFrom)
    rdTo$._oneToOne <- rownames(rdTo)

    ## Create the assay link
    hits <- .get_Hits(rdFrom, rdTo, "._oneToOne", "._oneToOne")
    al <- AssayLink(name = to,
                    from = from,
                    fcol = "._oneToOne",
                    hits = hits)

    ## Update the assay link in the QFeatures object
    .update_assay_links(object, al)
}
