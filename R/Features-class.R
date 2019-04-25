##' @title Quantitative MS Features
##'
##' @aliases Features Features-class class:Features
##'
##' @name Features
##'
##' @description
##'
##' Conceptually, a `Features` object holds a set of `matrix` (or
##' `array`) elements containing quantitative data. The number of
##' columns (samples) are always the same across the matrices, but the
##' number of rows (features) can vary. Each one of these matrices has
##' a set of feature annotation (encoded as `DataFrame` objects), that
##' have the same number of rows as the assay matrix their are
##' associated to, and an arbitrary number of columns (feature
##' variables). In addition, a `Features` object also uses a single
##' `DataFrame` to annotate the samples (columns) represented in all
##' the matrices. The largest assay matrix (the one with the highest
##' number of features) is considered the main assay, from which the
##' other ones are derived by aggregating/combining several rows into
##' a single one.
##'
##' A typical use case for such `Features` object if to represent
##' quantitative proteomics or metabolomics data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values.
##'
##' `Features` objects can be created with the `Features` constructor
##' or using the `[readFeatures()]` function to create an instance
##' from tabular data. The constructor, that is used to create objects
##' from their bare parts, that takes the following parameters:
##'
##' TODO
##'
##' @details
##'
##' More details come here ...
##'
##' @rdname Features-class
##'
##' @import S4Vectors
##' @importFrom SummarizedExperiment assays colData
##' @importFrom Biobase featureNames sampleNames
##'
##' @md
##'
##' @exportClass Features
##'
##' @author Laurent Gatto
##'
##' @examples
##'
##' ## Creating a Features object manually
##'
##' m1 <- matrix(1:40, ncol = 4)
##' m2 <- matrix(1:16, ncol = 4)
##' colnames(m1) <- colnames(m2) <- paste0("S", 1:4)
##' rownames(m1) <- letters[1:10]
##' rownames(m2) <- letters[1:4]
##'
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = letters[1:10])
##' df2 <- DataFrame(row.names = letters[1:4])
##'
##' x <- new("Features",
##'          assays = SimpleList(assay1 = m1, assay2 = m2),
##'          featureData = SimpleList(assay1 = df1, assay2 = df2),
##'          colData = DataFrame(row.names = colnames(m1)),
##'          metadata = list(paste("Generated on", Sys.Date())))
##' x
##'
##' ## Creating a Features object from a data.frame
##' data(hlpsms)
##' ft <- readFeatures(hlpsms, ecol = 1:10)
##' ft
##'
##' ## The assay isn't named yet, so let's name it 'psms', to clarify that the
##' ## data are PSM-level quantitations.
##' names(ft)
##' names(ft) <- "psms"
##' ft
NULL

setClass("Features",
    representation(
        assays = "SimpleList",
        featureData = "SimpleList",
        colData = "DataFrame",
        metadata = "list",
        version = "character"),
    prototype = prototype(version = "0.1"))

setValidity("Features", .valid_Features)

setMethod("show", "Features",
          function(object) {
              selectSome <- S4Vectors:::selectSome
              scat <- function(fmt, vals = character(), exdent = 2, ...) {
                  vals <- ifelse(nzchar(vals), vals, "''")
                  lbls <- paste(S4Vectors:::selectSome(vals), collapse = " ")
                  txt <- sprintf(fmt, length(vals), lbls)
                  cat(strwrap(txt, exdent = exdent, ...), sep="\n")
              }
              if (isEmpty(object)) .show_empty_Features(object)
              else .show_Features(object)
          })


setMethod("names", "Features", function(x) names(x@assays))
##' @exportMethod assays
setMethod("assays", "Features", function(x, ...) x@assays)
setMethod("dim", "Features", function(x) dim(x@assays[[main_assay(x)]]))
setMethod("length", "Features", function(x) length(x@assays))
setMethod("isEmpty", "Features", function(x) length(x) == 0)
##' exportMethod colData
setMethod("colData", "Features", function(x) x@colData)
setMethod("metadata", "Features", function(x, ...) x@metadata)
##' @importFrom Biobase featureData
setMethod("featureData", "Features", function(object) object@featureData)
##' @exportMethod featureNames
setMethod("featureNames", "Features",
          function(object) rownames(object@featureData[[main_assay(object)]]))
setGeneric("featureVariables", function(object) standardGeneric("featureVariables"))
##' @exportMethod featureVariables
setMethod("featureVariables", "Features",
          function(object) .featureVariables(object))

##' @exportMethod sampleNames
setMethod("sampleNames", "Features",
          function(object) rownames(object@colData))
setReplaceMethod("names", "Features",
    function(x, value) {
        names(x@assays) <- value
        names(x@featureData) <- value
        x
    })

setReplaceMethod("metadata", "Features",
    function(x, value) {
        if (!is.list(value))
            stop("replacement 'metadata' value must be a list")
        if (!length(value))
            names(value) <- NULL
        x@metadata <- value
        x
    })
