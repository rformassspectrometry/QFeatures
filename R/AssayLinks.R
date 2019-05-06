##' Links between Assays
##'
##' Links between assays within a [Features] object are handled by an
##' `AssayLinks` object. It is composed by a list of `AssayLink`
##' instances.
##'
##' The constructor are respectively `AssayLinks(..., .Names = NULL)`
##' and `AssayLink(name, from = NULL, fcol = NULL, hits = NULL)`. See
##' examples below.
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
##' AssayLinks(names = c("Assay1", "Assay2"))
NULL

##' @exportClass AssayLink
setClass("AssayLink",
         slots = c(name = "character",
                   from = "character",
                   fcol = "character",
                   hits = "Hits"))

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
##' @export
AssayLink <- function(name, from = NA_character_,
                      fcol = NA_character_,
                      hits = Hits()) {
    new("AssayLink", name = name[1],
        from = from, fcol = fcol,
        hits = hits)
}

##' @rdname AssayLinks
##' @param object An `AssayLink` object to show.
##' @export
setMethod("show", "AssayLink",
          function(object) cat("AssayLink for assay <", object@name, ">\n", sep = ""))

##' @exportClass AssayLinks
setClass("AssayLinks",
         contains = "SimpleList",
         prototype = prototype(
             elementType = "AssayLink"))

##' @param ... A set of `AssayLink` objects or a list thereof.
##' @param names A `character()` of `AssayLink` names. If provided,
##'     `...` are ignored, and `names` is used to create an
##'     `AssayLinks` object with `AssayLink` instances with names
##'     `names`.
##' @rdname AssayLinks
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


assayLink <- function(object, i)
    object@AssayLinks[[i]]

parentAssayLink <- function(x, i) {
    i2 <- assayLink(x, i)@from
    assayLink(x, i2)
}


