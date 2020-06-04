##' @title Links between Assays
##'
##' @description
##'
##' Links between assays within a [Features] object are handled by an
##' `AssayLinks` object. It is composed by a list of `AssayLink`
##' instances.
##'
##' @section Constructors:
##'
##' Object can be created with the `AssayLink()` and `AssayLinks()` constructors.
##'
##' @section Methods and functions:
##'
##' - `assayLink(x, i)` accesses the AssayLink at position `i` or with
##'    name `i` in the [Features] object `x`.
##'
##' - `parentAssayLinks(x, i, recursive = FALSE)` accesses the
##'   parent(s) `AssayLinks` or assay with index or name `i`.
##' 
##' @section Creating links between assays:
##' 
##' - `addAssayLink` takes any two assays from the [Features] object and 
##'   creates a link given a matching feature variable in each assay's 
##'   `rowData`.
##' - `addAssayLinkOneToOne` also links two assays from the [Features] 
##'   object, but the assays must have the same size and contain the same 
##'   rownames although a different ordering is allowed. The matching is 
##'   performed based on the row names of the assays. 
##' - `addAssayLinkMultiParent` links a child assay to multiple parent assays in 
##'   a [Features] object. The matching is performed based on the supplied 
##'   variables names contained in the assay's `rowData`.
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
##' ## create a Features object with 2 (identical) assays
##' ## see also '?Features'
##' se <- SummarizedExperiment(matrix(runif(20), ncol = 2, 
##'                                   dimnames = list(LETTERS[1:10], letters[1:2])),
##'                            rowData = DataFrame(ID = 1:10))
##' ft <- Features(list(assay1 = se, assay2 = se))
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
##' ## An assay can also be linked to multiple parent assays
##' ## create a Features object with 2 parent assays and 1 child assay
##' ft <- Features(list(parent1 = se[1:6, ], parent2 = se[4:10, ], child = se))
##' ft <- addAssayLinkMultiParent(ft, from = c("parent1", "parent2"),
##'                               to = "child", varsFrom = c("ID", "ID"),
##'                               varTo = "ID")
##' assayLink(ft, "child")
##'
NULL

setClassUnion("ListOrHits", c("Hits", "List"))

##' @exportClass AssayLink
setClass("AssayLink",
         slots = c(name = "character",
                   from = "character",
                   fcol = "character",
                   hits = "ListOrHits"))

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
                                   paste(sapply(object@hits, length), collapse = ","),
                                   length(object@hits)), "]\n", sep = "")
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
    names(args) <- sapply(args, slot, "name")    
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


##' @param x An instance of class [Features].
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


##' @param j ignored.
##' @param drop ignored.
##' @rdname AssayLinks
setMethod("[", c("AssayLink", "character"),
          function(x, i, j, ..., drop = TRUE) {
              k <- which(elementMetadata(x@hits)$names_to %in% i)
              x@hits <-  x@hits[k, ]
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

## The function takes the rowData of an assays to link from, the rowData of an
## assay to link to, the corresponding assay names and the feature variable 
## names in both rowData that link the 2 assays together. The function returns 
## an AssayLink object that links the two assays given the relationship between 
## the corresponding feature variables.
.create_assay_link <- function(rdFrom, 
                               rdTo,
                               from,
                               to,
                               varFrom, 
                               varTo) {
    if (identical(from, to))
        stop("Creating an AssayLink between an assay and itself is not allowed.")
    ## Get the shared feature variable 
    matchFrom <- unlist(rdFrom[varFrom], use.names = FALSE)
    matchTo <- unlist(rdTo[varTo], use.names = FALSE)
    ## Find hits
    hits <- findMatches(matchFrom, matchTo)
    if (length(c(hits@from, hits@to)) == 0) 
        stop(paste0("No match found between field '", varFrom, "' (in '", 
                    from, "') and filed '", varTo, "' (in '", to, "')."))
    ## Add the row names corresponding to the hits
    elementMetadata(hits)$names_from <- rownames(rdFrom)[hits@from]
    elementMetadata(hits)$names_to <- rownames(rdTo)[hits@to]
    ## Create the new link
    AssayLink(name = to,
              from = from,
              fcol = varFrom,
              hits = hits)
}


## The function takes a Features object, the assay name*s* to link from, the 
## assay name to link to and the feature variable names contained in the rowData 
## of the assay*s* to link from and the assay to link to. The function returns 
## an AssayLink object that links the parent assay*s* to the child assay given 
## the relationship between the corresponding feature variables.
.create_assay_link_multi <- function(object, 
                                     from, 
                                     to, 
                                     varsFrom, 
                                     varTo) {
    if (missing(varsFrom) | missing(varTo)) {
        ## Create the list of hits between the child assay and the parent assays
        ## based on the rownames 
        rdTo <- cbind(rowData(object[[to]]), RowNames = rownames(object[[to]]))
        hitsList <- lapply(seq_along(from), function(ii) {
            rdFrom <- cbind(rowData(object[[from[ii]]]), 
                            RowNames = rownames(object[[from[[ii]]]]))
            .create_assay_link(rdFrom = rdFrom,
                               rdTo = rdTo,
                               from[ii], to,
                               "RowNames", "RowNames")@hits
        })
        varsFrom <- "RowNames"
    } else if (length(from) == length(varsFrom)) {
        ## Create the list of hits between the child assay and the parent assays
        ## based on the supplied varsFrom and varTo
        hitsList <- lapply(seq_along(from), function(ii) {
            .create_assay_link(rdFrom = rowData(object[[from[ii]]]),
                               rdTo = rowData(object[[to]]),
                               from[ii], to,
                               varsFrom[[ii]], varTo)@hits
        })
    } else
        stop("'from' and 'varsFrom' must have same length.")
    
    names(hitsList) <- from
    ## Return a multi-parent AssayLink 
    AssayLink(name = to, 
              from = from, 
              fcol = varsFrom,
              hits = List(hitsList))
}



## Function that updates the Features object's AssayLinks with the provided 
## AssayLink object. 
.update_assay_links <- function (object, al) {
    if (!inherits(al ,"AssayLink")) stop("'al' must be an AssayLink object.")
    ## Check the child indexing on rownames
    if (!all(elementMetadata(al@hits)$names_to %in% rownames(object[[al@name]])))
        stop("Invalid AssayLink. The AssayLink metadata 'names_to' does not match the rownames.")
    ## Check the parent indexing on rownames 
    
    if (!all(elementMetadata(al@hits)$names_from %in% rownames(object[[al@from]])))
        stop("Invalid AssayLink. The AssayLink metadata 'names_from' does not match the rownames.")
    
    ## Overwrite the AssayLink
    object@assayLinks@listData[[al@name]] <- al
    
    if (validObject(object))
        return(object)
}

## Same as '.update_assay_links' but when the AssayLink contains multiple 
## parents 
.update_assay_links_multi_parents <- function (object, al) {
    if (!inherits(al ,"AssayLink")) stop("'al' must be an AssayLink object.")
    ## Check the child indexing on rownames
    isCorrectToLink <- sapply(al@hits, function(l) all(elementMetadata(l)$names_to %in% rownames(object[[al@name]])))
    if (any(!isCorrectToLink)) 
        stop("Invalid AssayLink. At least one of the 'hits' metadata 'names_to' does not match the rownames.")
    ## Check the parent indexing on rownames 
    isCorrectFromLink <- sapply(seq_along(al@hits), 
                              function(i) all(elementMetadata(al@hits[[i]])$names_from %in% rownames(object[[al@from[i]]])))
    if (any(!isCorrectFromLink))
        stop("Invalid AssayLink. The AssayLink metadata 'names_from' does not match the rownames.")
    
    ## Overwrite the AssayLink
    object@assayLinks@listData[[al@name]] <- al
    
    if (validObject(object))
        return(object)
}

##' @rdname AssayLinks
##'
##' @param from A character(1) or integer(1) indicating which assay to link from
##'     in `object`
##' @param to A character(1) or integer(1) indicating which assay to link to in
##'     `object`
##' @param varFrom A character (1) indicating the feature variable to use to 
##'     match the `from` assay to the `to` assay. 
##' @param varTo A character (1) indicating the feature variable to use to 
##'     match the `to` assay to the `from` assay.
##'
##' @export
addAssayLink <- function(object, 
                         from, 
                         to,
                         varFrom, 
                         varTo) {
    if (is.numeric(from)) from <- names(object)[[from]]
    if (is.numeric(to)) to <- names(object)[[to]]
    ## Create the assay link
    al <- .create_assay_link(rdFrom = rowData(object[[from]]),
                             rdTo = rowData(object[[to]]),
                             from, to,
                             varFrom, varTo)
    ## Update the assay link in the Features object
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
    rdFrom$OneToOne <- rownames(rdFrom)
    rdTo$OneToOne <- rownames(rdTo)
    ## Create the assay link
    al <- .create_assay_link(rdFrom, rdTo,
                             from, to, 
                             "OneToOne", "OneToOne")
    ## Update the assay link in the Features object
    .update_assay_links(object, al)
}

##' @rdname AssayLinks
##' 
##' @param varsFrom A `character()` indicating the feature variable**s** to use 
##'     to match the `from` assay**s** to the `to` assay. `varsFrom` is assumed 
##'     to be ordered as `from`, and both arguments must have same length. If 
##'     `length(varsFrom) == length(from) == 1`, then `addAssayLink` is called
##'     instead. 
##' 
##' @export
addAssayLinkMultiParent <- function(object, 
                                    from, 
                                    to, 
                                    varsFrom,
                                    varTo) {
    if (length(from) == 1) {
        warning("Only 1 parent supplied, calling 'addAssay'")
        return(addAssayLink(object, from, to, varsFrom, varTo))
    }
    if (is.numeric(from)) from <- names(object)[[from]]
    if (is.numeric(to)) to <- names(object)[[to]]
    
    al <- .create_assay_link_multi(object, from, to, varsFrom, varTo)                  
    
    ## Update the assay link in the Features object
    .update_assay_links_multi_parents(object, al)
}
