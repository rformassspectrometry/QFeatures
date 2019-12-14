.zeroIsNA <- function(object) {
    sel <- assay(object) == 0
    assay(object)[sel] <- NA
    object
}

.nNA <- function(object) {
    nNA <- sum(is.na(assay(object)))/prod(dim(object))
    nNA_rows <- rowSums(is.na(assay(object)))
    nNA_cols <- colSums(is.na(assay(object)))
    return(nNA = nNA, nNArows = nNA_rows, nNAcols = nNA_cols)
}


##' @title Managing missing data
##'
##' @description
##' 
##' This manual page describes the handling of missing values in
##' `Features` objects. The current functions are available:
##'
##' - `zeroIsNA(object)`: this function replaces all 0 in `object` by
##'    `NA`. This is often necessary when third-party software assume
##'    that features that weren't quantified should be assigned an
##'    intensity of 0.
##'
##' 
##' @param object An object of class `Features` or `SummarizedExperiment`.
##'
##' @param i The index or name of the assay to be processed.
##'
##' @return An instance of the same class as `object`.
##'
##' @aliases zeroIsNA zeroIsNA,SummarizedExperiment,missing-method zeroIsNA,Features,missing-method zeroIsNA,Features,numeric-method zeroIsNA,Features,integer-method zeroIsNA,Features,character-method
##'
##' @name missing-data
##' 
##' @rdname missing-data
NULL


##' @exportMethod zeroIsNA
setMethod("zeroIsNA", c("SummarizedExperiment", "missing"),
          function(object, i) .zeroIsNA(object))

setMethod("zeroIsNA", c("Features", "missing"),
          function(object, i) 
              for (i in seq_len(length(object)))
                  object[[i]] <- zeroIsNA(object[[i]]))

setMethod("zeroIsNA", c("Features", "numeric"),
          function(object, i) 
                  object[[i]] <- zeroIsNA(object[[i]]))

setMethod("zeroIsNA", c("Features", "integer"),
          function(object, i) 
                  object[[i]] <- zeroIsNA(object[[i]]))

setMethod("zeroIsNA", c("Features", "character"),
          function(object, i) 
                  object[[i]] <- zeroIsNA(object[[i]]))



