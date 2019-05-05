##' Quantitative MS Features
##' 
##' Conceptually, a `Features` object holds a set of *assays*, each
##' composed of a `matrix` (or `array`) containing quantitative data
##' and row annotations (meta-data).  The number and the names
##' (nature) of columns (samples) are always the same across the
##' assays, but the number of rows (features) can vary. The assays are
##' typically defined as `MSnSet` or `SummarizedExperiment`
##' objects. In addition, a `Features` object also uses a single
##' `DataFrame` to annotate the samples (columns) represented in all
##' the matrices.
##'
##' The `Features` class extends the
##' [MultiAssayExperiment::MatchedAssayExperiment] and inherits all
##' the functionality of the
##' [MultiAssayExperiment::MultiAssayExperiment] class.
##'
##' A typical use case for such `Features` object is to represent
##' quantitative proteomics (or metabolomics) data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values.
##'
##' The largest assay (the one with the highest number of features) is
##' considered the main assay, from which the other ones are derived
##' by aggregating/combining several rows into a single one.
##'
##' The recommended way to create `Features` objects is the use the
##' `readFeatures()` function, that creates an instance from tabular
##' data. The `Features` constructor can be used to create objects
##' from their bare parts.  It is the user's responsability to make
##' sure that these match the class validity requirements.
##'
##'
##' @section Accessors:
##'
##' See [MultiAssayExperiment()].
##'
##' @section Working with Features:
##'
##' - `addAssay(object, x, name, assayLinks)`: Adds a new assay `x` to
##'   the `Features` instance `object`. 
##'
##' @seealso The [readFeatures()] constructor.
##'
##' @import MultiAssayExperiment
##'
##' @name Features
##'
##' @rdname Features-class
##'
##' @aliases Features Features-class class:Features addAssay
##'
##' @md
##'
##' @exportClass Features
##'
##' @author Laurent Gatto
##'
##' @examples
##' ## ------------------------
##' ## An empty Features object
##' ## ------------------------
##' 
##' Features()
##'
##' ## -----------------------------------
##' ## Creating a Features object manually
##' ## -----------------------------------
##'
##' ## two assays (matrices) with matching column names
##' m1 <- matrix(1:40, ncol = 4)
##' m2 <- matrix(1:16, ncol = 4)
##' sample_names <- paste0("S", 1:4)
##' colnames(m1) <- colnames(m2) <- sample_names
##' rownames(m1) <- letters[1:10]
##' rownames(m2) <- letters[1:4]
##'
##' ## two corresponding feature metadata with appropriate row names
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = rownames(m1))
##' df2 <- DataFrame(row.names = rownames(m2))
##'
##' (se1 <- SummarizedExperiment(m1, df1))
##' (se2 <- SummarizedExperiment(m2, df2))
##'
##' ## Sample annotation (colData)
##' cd <- DataFrame(Var1 = rnorm(4),
##'                 Var2 = LETTERS[1:4],
##'                 row.names = sample_names)
##'
##' el <- list(assay1 = se1, assay2 = se2)
##' fts1 <- Features(el, colData = cd)
##' fts1
##'
##' ## Add an assay
##' fts1 <- addAssay(fts1, se1[1:2, ], name = "se3")
##'
##' ## -------------------------------------------------
##' ## Creating a Features object from a data.frame (see
##' ## ?readFeatures) for details
##' ## -------------------------------------------------
##' 
##' data(hlpsms)
##' fts2 <- readFeatures(hlpsms, ecol = 1:10, name = "psms")
##' fts2
##'
##' fts2[[1]]
##' fts2[["psms"]]
NULL

setClass("Features",
         contains = "MatchedAssayExperiment",
         slots = c(version = "character",
                   links = "DataFrame"),
         prototype = prototype(
             version = "0.1",
             links = EmptyAssayLinks()))


setMethod("show", "Features",
          function(object) {
              if (isEmpty(object)) {
                  cat(sprintf("A empty instance of class %s", class(object)), "\n")
                  return(NULL)
              }
              cat(sprintf("A instance of class %s", class(object)), "containing")
              el <- experiments(object)
              o_class <- class(el)
              elem_cl <- vapply(el, class, character(1L))
              o_len <- length(el)
              o_names <- names(el)
              featdim <- vapply(el, FUN = function(obj) {
                  dim(obj)[1]
              }, FUN.VALUE = integer(1L))
              sampdim <- vapply(el, FUN = function(obj) {
                  dim(obj)[2]
              }, FUN.VALUE = integer(1L))
              cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                          seq(o_len), o_names, elem_cl, featdim, sampdim), "\n")              
          })


