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
##' al1 <- AssayLink(name = "assay1")
##' al1
##' 
##' AssayLinks(al1)
##'
##' al2 <- AssayLinks(names = c("Assay1", "Assay2"))
##' al2
NULL

##' @exportClass AssayLink
setClass("AssayLink",
         slots = c(name = "character",
                   from = "character",
                   fcol = "character",
                   hits = "Hits"))

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
              cat("AssayLink for assay <", object@name, ">\n", sep = "")
              cat("[from:", object@from, "|fcol:", object@fcol, "|hits:", length(object@hits),"]\n", sep = "")
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
    if (is.na(this@from))
        return(AssayLinks(this))
    ans <- list()
    while (!is.na(this@from)) {
        ans <- append(ans, this)
        this <- assayLink(x, this@from)
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

## The function takes a Features object, two assay names and two feature 
## variable names and creates an AssayLink object that links the two assays.
.createAssayLink <- function (object, ## Features object
                             from, 
                             to,
                             varFrom, 
                             varTo) {
    ## Get the shared feature variable 
    rowDatFrom <- unlist(rowData(object[[from]])[varFrom], use.names = FALSE)
    rowDatTo <- unlist(rowData(object[[to]])[varTo], use.names = FALSE)
    ## Find hits
    hits <- findMatches(rowDatFrom, rowDatTo)
    if(length(c(hits@from, hits@to)) == 0) 
        stop(paste0("No match found between field '", varFrom, "' (in '", from, 
                    "') and filed '", varTo, "' (in '", to, "')."))
    ## Add the row names corresponding to the hits
    elementMetadata(hits)$names_from <- rownames(object[[from]])[hits@from]
    elementMetadata(hits)$names_to <- rownames(object[[to]])[hits@to]
    ## Create the new link
    al <- AssayLink(name = to,
                    from = from,
                    fcol = varFrom,
                    hits = hits)
}


#' @rdname AssayLinks
#'
#' @param from A character(1) or integer(1) indicating which assay to link from
#'     in `object`
#' @param to A character(1) or integer(1) indicating which assay to link to in
#'     `object`
#' @param varFrom A character (1) indicating the feature variable to use to 
#'     match the `from` assay to the `to` assay. 
#' @param varTo A character (1) indicating the feature variable to use to 
#'     match the `to` assay to the `from` assay. If missing, `varTo` is the 
#'     same as `varFrom`.
#'
#' @export
createAssayLink <- function(object, 
                            from, 
                            to,
                            varFrom, 
                            varTo){
    if (missing(varTo)) varTo <- varFrom
    if (is.numeric(from)) from <- names(object)[[from]]
    if (is.numeric(to)) to <- names(object)[[to]]
    al <- .createAssayLink(object, from, to, varFrom, varTo)
    ## Append assay link to the existing links
    object@assayLinks@listData[[to]] <- al
    stopifnot(validObject(object))
    object
}

createOneToOneAssayLink <- function(object, 
                                    from, 
                                    to, 
                                    fcol){
    ## Get number of rows for each assay to link
    N <- unique(dims(object)[1, c(from, to)])
    if (length(N) != 1) 
        stop("The 'from' and 'to' assays must have the same number of rows.")
    if (missing(fcol)) {
        ## Create a dummy variable that makes a 1-1 link between asssays
        rowData(object[[from]])$oneToOneID <- 1:N
        rowData(object[[to]])$oneToOneID <- 1:N
        fcol <- "oneToOneID"
    }
    ## Create the 1-1 link
    createAssayLink(object, from, to, fcol)
}

