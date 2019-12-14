.zeroIsNA <- function(x) {
    sel <- assay(x) == 0
    assay(x)[sel] <- NA
    x
}

.nNA <- function(x) {
    nNA <- sum(is.na(assay(x)))/prod(dim(x))
    nNA_rows <- table(rowSums(is.na(assay(x))))
    nNA_cols <- colSums(is.na(assay(x)))
    return(list(nNA = nNA, nNArows = nNA_rows, nNAcols = nNA_cols))
}


##' @title Managing missing data
##'
##' @description
##' 
##' This manual page describes the handling of missing values in
##' `Features` objects. In the following functions, if `object` is of
##' class `Feature`, and optional assay index of name `i` can be
##' specified to define the assay (by name of index) on which to
##' operate.
##'
##' The following functions are currently available:
##'
##' - `zeroIsNA(object, i)` replaces all 0 in `object` by `NA`. This
##'    is often necessary when third-party software assume that
##'    features that weren't quantified should be assigned an
##'    intensity of 0. 
##'
##' - `nNA(object)` return a list of missing value summaries. The
##'   first element `nNA` gives the percentage of missing values; the
##'   second element `nNArows` provides a table of the number of
##'   missing values for the features (rows) of the assay(s); the
##'   third element `nNAcols` provides the number of missing values in
##'   each sample of the assay(s).
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
##' @aliases nNA nNA,SummarizedExperiment,missing-method nNA,Features,missing-method nNA,Features,numeric-method nNA,Features,integer-method nNA,Features,character-method
##'
##' @name missing-data
##' 
##' @rdname missing-data
NULL


##' @exportMethod zeroIsNA
setMethod("zeroIsNA", c("SummarizedExperiment", "missing"),
          function(object, i) .zeroIsNA(object))

setMethod("zeroIsNA", c("Features", "missing"),
          function(object, i) {
              for (i in seq_len(length(object)))
                  object[[i]] <- zeroIsNA(object[[i]])
              object
          })

setMethod("zeroIsNA", c("Features", "integer"),
          function(object, i) {
              object[[i]] <- zeroIsNA(object[[i]])
              object
          })

setMethod("zeroIsNA", c("Features", "numeric"),
          function(object, i) zeroIsNA(object, as.integer(i)))

setMethod("zeroIsNA", c("Features", "character"),
          function(object, i) {
              object[[i]] <- zeroIsNA(object[[i]])
              object
          })


##' @exportMethod nNA
setMethod("nNA", c("SummarizedExperiment", "missing"),
          function(object, i) .nNA(object))

setMethod("nNA", c("Features", "integer"),
          function(object, i) .nNA(object[[i]]))

setMethod("nNA", c("Features", "numeric"),
          function(object, i) .nNA(object[[as.integer(i)]]))

setMethod("nNA", c("Features", "character"),
          function(object, i) .nNA(object[[i]]))

setMethod("nNA", c("Features", "missing"),          
          function(object, i) {
              if (length(object) == 1)
                  return(nNA(object, 1))
              res <- lapply(seq_len(length(object)),
                            function(i) .nNA(object[[i]]))
              names(res) <- names(object)
              ans <- vector("list", length = 3)
              names(ans) <- c("nNA", "nNArows", "nNAcols")
              ans[[1]] <- sapply(res, "[[", 1)              
              ans[[3]] <- t(sapply(res, "[[", 3))
              ans2 <- matrix(0,
                             ncol = 1 + nrow(colData(object)),
                             nrow = ncol(colData(object)))
              rownames(ans2) <- names(object)
              colnames(ans2) <- 0:nrow(colData(object))
              for (i in seq_len(length(res))) {
                  x <- res[[i]]$nNArows
                  ans2[i, names(x)] <- x
              }
              ans[[2]] <- ans2
              ans
          })



