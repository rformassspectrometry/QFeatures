##' Links between Assays
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
##' @aliases AssayLinks AssayLink AssayLink-class AssayLinks-class class:AssayLinks class:AssayLink show,AssayLink-method
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
##'     `combineFeatures`). `NA_character_`, if not applicable.
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
    object@assayLinks[[i]]


##' @param x An instance of class [Features].
##' 
##' @param i The index or name of the assay whose `AssayLink` or
##'     parent's are to be returned.
##' 
##' @param recursive If `TRUE`, returns all parents. Default is
##'     `FALSE`.
##' 
##' @return `parentAssayLink` returns an instance of class
##'     `AssayLinks` or `NULL`, if no parent available. 
##' 
##' @export
##' @md
##' @rdname AssayLinks
parentAssayLinks <- function(x, i, recursive = FALSE) {
    i2 <- assayLink(x, i)@from
    if (is.na(i2))
        return(NULL)
    par0 <- par <- assayLink(x, i2)
    if (!recursive) par0 <- NULL
    while (!is.null(par0)) {
        par0 <- assayLink(x, par0@from)
        par <- append(par, par0)
    }
    return(AssayLinks(par))
}
