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
##' `createAssayLink` and `createAssayLinkOneToOne` link two assays in a 
##' Features object. 
##' 
##' - `createAssayLink` takes any two assays from the `Features` object and 
##'   creates a link given a matching feature variable in each assay's 
##'   `rowData`.
##' - `createAssayLinkOneToOne` also links two assays from the `Features` 
##'   object, but the assays must have the same size and a one-to-one 
##'   relationship between the assays is created. A common feature variable can 
##'   be supplied, otherwise the rows of the assays are assumed to be matched.
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

## The function takes two SummarizedExperiments objects, the corresponding assay 
## names and the feature variable and creates an AssayLink object that links the 
## two assays given the .
.createAssayLink <- function (seFrom, 
                              seTo,
                              nameFrom,
                              nameTo,
                              varFrom, 
                              varTo) {
    if(identical(nameFrom, nameTo))
        stop("Creating an AssayLink between an assay and itself is not allowed.")
    ## Get the shared feature variable 
    rowDatFrom <- unlist(rowData(seFrom)[varFrom], use.names = FALSE)
    rowDatTo <- unlist(rowData(seTo)[varTo], use.names = FALSE)
    ## Find hits
    hits <- findMatches(rowDatFrom, rowDatTo)
    if(length(c(hits@from, hits@to)) == 0) 
        stop(paste0("No match found between field '", varFrom, "' (in '", 
                    nameFrom, "') and filed '", varTo, "' (in '", nameTo, "')."))
    ## Add the row names corresponding to the hits
    elementMetadata(hits)$names_from <- rownames(seFrom)[hits@from]
    elementMetadata(hits)$names_to <- rownames(seTo)[hits@to]
    ## Create the new link
    al <- AssayLink(name = nameTo,
                    from = nameFrom,
                    fcol = varFrom,
                    hits = hits)
}

## The function takes two SummarizedExperiments objects and creates an AssayLink 
## object that links the two assays with a one-to-one relation.
.createAssayLinkOneToOne <- function (seFrom,
                                      seTo,
                                      nameFrom,
                                      nameTo,
                                      varCommon) {
    ## Get number of rows for each assay to link
    N <- unique(c(nrow(seFrom), nrow(seTo)))
    if (length(N) != 1) 
        stop("The 'from' and 'to' assays must have the same number of rows.")
    if (missing(varCommon)) {
        ## Create a dummy variable that makes a 1-1 link between asssays
        rowData(seFrom)$oneToOneID <- 1:N
        rowData(seTo)$oneToOneID <- 1:N
        varCommon <- "oneToOneID"
    }
    ## Create the 1-1 link
    .createAssayLink(seFrom, seTo, 
                     nameFrom, nameTo, 
                     varCommon, varCommon)
}

## Function that adds the provided AssayLink object to the provided Features 
## object
## TODO adapt this when allowing an assay to have several parents 
## TODO would it be useful for the end user to have this exported 
.addAssayLink <- function (object, al) {
    if(!inherits(al ,"AssayLink")) stop("'al' must be an AssayLink object.")
    object@assayLinks@listData[[al@name]] <- al
    stopifnot(validObject(object))
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
##'     match the `to` assay to the `from` assay. If missing, `varTo` is the 
##'     same as `varFrom`.
##'
##' @export
createAssayLink <- function (object, 
                             from, 
                             to,
                             varFrom, 
                             varTo) {
    if (missing(varTo)) varTo <- varFrom
    if (is.numeric(from)) from <- names(object)[[from]]
    if (is.numeric(to)) to <- names(object)[[to]]
    al <- .createAssayLink(seFrom = object[[from]],
                           seTo = object[[to]],
                           nameFrom = from,
                           nameTo = to,
                           varFrom, varTo)
    .addAssayLink(object, al)
}

##' @rdname AssayLinks
##' 
##' @param varCommon The name of the common feature variable in `from` and `to`. 
##'     If missing, row indexing will be used as matching variable.
##' 
##' @export
createAssayLinkOneToOne <- function (object, 
                                     from, 
                                     to, 
                                     varCommon) {
    if (is.numeric(from)) from <- names(object)[[from]]
    if (is.numeric(to)) to <- names(object)[[to]]
    al <- .createAssayLinkOneToOne(seFrom = object[[from]],
                                   seTo = object[[to]],
                                   nameFrom = from,
                                   nameTo = to,
                                   varCommon = varCommon)
    .addAssayLink(object, al)
}

