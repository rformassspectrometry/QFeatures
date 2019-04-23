setClass("Features",
    representation(
        assays = "SimpleList",
        fData = "DataFrame",
        colData = "DataFrame",
        metadata = "list"
    ))

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


setMethod("assayNames", "Features", function(x, ...) names(x@assays))

setMethod("names", "Features", function(x) names(x@assays))

setMethod("assays", "Features", function(x, ...) x@assays)


setMethod("dim", "Features", function(x) dim(x@assays[[.main_assay(x)]]))

setMethod("length", "Features", function(x) length(x@assays))

setMethod("isEmpty", "Features", function(x) length(x) == 0)

setMethod("colData", "Features", function(x) x@colData)

setMethod("metadata", "Features", function(x, ...) x@metadata)

setMethod("featureNames", "Features",
          function(object) rownames(object@fData))

setMethod("sampleNames", "Features",
          function(object) rownames(object@colData))

setMethod("dimnames", "Features",
    function(x) {
    list(featureNames(x), rownames(colData(x)))
})

setReplaceMethod("names", "Features",
    function(x, value) {
        names(x@assays) <- value
        x
    })

setReplaceMethod("assayNames", "Features",
    function(x, value) {
        names(x@assays) <- value
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

Features <- function(assays = SimpleList(),
                     fData = DataFrame(),
                     colData = DataFrame()) {
    new("Features",
        assays = assays,
        fData = fData,
        colData = colData)
}
