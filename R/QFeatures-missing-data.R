###---- Internal functions ----####

.zeroIsNA <- function(x) {
    sel <- assay(x) == 0
    assay(x)[sel] <- NA
    x
}

.infIsNA <- function(x) {
    sel <- is.infinite(assay(x))
    assay(x)[sel] <- NA
    x
}

## Internal function that compute the number and proportion of missing
## data from a SummarizedExperiment object
.nNAByAssay <- function(object) {
    x <- assay(object)
    if (nrow(x) == 0 || ncol(x) == 0) {
        nNA <- NA_integer_
        pNA <- NA_real_
    } else {
        nNA <- sum(is.na(x))
        pNA <- nNA / length(x)
    }
    DataFrame(nNA = unname(nNA),
              pNA = unname(pNA))
}

## Internal function that compute the number and proportion of missing
## data for a given margin (features = 1 and sample = 2) from a
## SummarizedExperiment object
.nNAByMargin <- function(object, MARGIN = 1) {
    x <- assay(object)
    if (nrow(x) == 0 || ncol(x) == 0) {
        nNA <- NA_integer_
        pNA <- NA_real_
        names(pNA) <- ""
    } else {
        nNA <- apply(is.na(x), MARGIN, sum)
        n <- ifelse(MARGIN == 1, ncol(x), nrow(x))
        pNA <- nNA / n
    }
    DataFrame(name = names(pNA),
              nNA = unname(nNA),
              pNA = unname(pNA),
              row.names = names(pNA))
}

## Internal function that compute the number and proportion of missing
## data for a SummarizedExperiment object
.nNA <- function(x, addToObject) {
    nNA <- .nNAByAssay(x)
    nNA_rows <- .nNAByMargin(x, 1)
    nNA_cols <- .nNAByMargin(x, 2)
    if (addToObject) {
        colData(x) <- cbind(colData(x), nNA_cols[, c("nNA", "pNA")])
        rowData(x) <- cbind(rowData(x), nNA_rows[, c("nNA", "pNA")])
        return(x)
    }
    list(nNA = nNA, nNArows = nNA_rows, nNAcols = nNA_cols)
}

## Internal function that compute the number and proportion of missing
## data for a QFeatures object
.nNAi <- function(object, i, addToObject) {
    i <- .normIndex(object, i)

    nNAcol <- do.call(rbind, lapply(i, function(ii)
        cbind(assay = ii, .nNAByMargin(object[[ii]], 2))
    ))

    if (!addToObject) {
        nNAassay <- do.call(rbind, lapply(i, function(ii)
            cbind(assay = ii, .nNAByAssay(object[[ii]]))
        ))

        nNArow <- do.call(rbind, lapply(i, function(ii)
            cbind(assay = ii, .nNAByMargin(object[[ii]], 1))
        ))

        return(list(
            nNA = nNAassay,
            nNArows = nNArow,
            nNAcols = nNAcol
        ))
    }

    for (ii in i) {
        se <- object[[ii]]
        rowStats <- .nNAByMargin(se, 1)[, c("nNA", "pNA")]
        rowData(se) <- cbind(rowData(se), rowStats)
        object[[ii]] <- se
    }

    cd <- colData(object)
    key <- rownames(cd)

    idx <- match(key, nNAcol$name)

    cd$nNA <- nNAcol$nNA[idx]
    cd$pNA <- nNAcol$pNA[idx]

    colData(object) <- cd
    object
}


.row_for_filterNA <- function(x, pNA = 0L) {
    if (!is.matrix(x))
        stop(sQuote("x"), " must be a matrix.")
    if (!is.numeric(pNA))
        stop(sQuote("pNA"), " must be numeric.")
    if (length(pNA) > 1)
        stop(sQuote("pNA"), " must be of length one.")
    if (pNA < 0 | pNA > 1)
        stop(sQuote("pNA"), " must be between 0 and 1.")
    k <- rowSums(is.na(x)) / ncol(x)
    k <= pNA
}


###---- Documentation ----####


##' @title Managing missing data
##'
##' @description
##'
##' This manual page describes the handling of missing values in
##' [QFeatures] objects. In the following functions, if `object` is of
##' class `QFeatures`, an optional assay index or name `i` can be
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
##' - `infIsNA(object, i)` replaces all infinite values in `object` by
##'    `NA`. This is necessary when third-party software divide
##'    expression data by zero values, for instance during custom
##'    normalization.
##'
##' - `nNA(object, i, addToObject = FALSE)` computes summaries of missing
##'   values. When `addToObject = FALSE`, it returns a list with three
##'   elements: (1) `nNA`, a `DataFrame` reporting the number and
##'   proportion of missing values for the whole set; (2) `nNArows`,
##'   a `DataFrame` with the number and proportion of missing values for
##'   the features (rows) of the set(s); and (3) `nNAcols`, a
##'   `DataFrame` with the number and proportion of missing values for
##'   each sample (columns) of the set(s). When `object` has class
##'   `QFeatures`, an additional column indicating the assay is included
##'   in each `DataFrame`. When `addToObject = TRUE`, no list is returned;
##'   instead, two columns, `nNA` (number of missing values) and `pNA`
##'   (proportion of missing values),  are added to the `rowData` and
##'   `colData` of the set(s).
##'
##' - `filterNA(object, pNA, i)` removes features (rows) that contain
##'   a proportion of more missing values of `pNA` or higher.
##'
##' See the *Processing* vignette for examples.
##'
##' @param object An object of class `QFeatures` or `SummarizedExperiment`.
##'
##' @param pNA `numeric(1)` providing the maximum proportion of
##'     missing values per feature (row) that is acceptable. Feature
##'     with higher proportions are removed. If 0 (default), features
##'     that contain any number of `NA` values are dropped.
##'
##' @param i One or more indices or names of the set(s) to be processed.
##'
##' @param addToObject `logical(1)` indicating if the nNA should
##'     return a list of `DataFrame` that contain the metrics or 
##'     a `QFeatures` object that contains the metrics in the 
##'     `rowData` and `colData` of the specified set(s).
##'
##' @return An instance of the same class as `object`.
##'
##' @aliases zeroIsNA zeroIsNA,SummarizedExperiment,missing-method zeroIsNA,QFeatures,missing-method zeroIsNA,QFeatures,numeric-method zeroIsNA,QFeatures,integer-method zeroIsNA,QFeatures,character-method
##'
##' @aliases infIsNA infIsNA,SummarizedExperiment,missing-method infIsNA,QFeatures,missing-method infIsNA,QFeatures,numeric-method infIsNA,QFeatures,integer-method infIsNA,QFeatures,character-method
##'
##' @aliases nNA nNA,SummarizedExperiment,missing-method nNA,QFeatures,missing-method nNA,QFeatures,numeric-method nNA,QFeatures,integer-method nNA,QFeatures,character-method
##'
##' @aliases filterNA filterNA,SummarizedExperiment-method filterNA,QFeatures-method
##'
##' @name missing-data
##'
##' @rdname QFeatures-missing-data
##'
##' @seealso The `impute()` for `QFeatures` instances.
##'
##' @examples
##' data(ft_na)
##'
##' ## Summary if missing values
##' nNA(ft_na, 1)
##'
##' ## Insert NA metrics into the QFeatures' rowData and colData
##'
##' ft_na <- nNA(ft_na, 1, addToObject = TRUE)
##'
##' ## Remove rows with missing values
##' assay(filterNA(ft_na, i = 1))
##'
##' ## Replace NAs by zero and back
##' ft_na <- impute(ft_na, method = "zero", i = 1, name = "imputedSet")
##' assay(ft_na[["imputedSet"]])
##' ## Replace zero by NA in the newly created set
##' ft_na <- zeroIsNA(ft_na, 2)
##' assay(ft_na[["imputedSet"]])
NULL


####---- zeroIsNA ----####


##' @exportMethod zeroIsNA
##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("SummarizedExperiment", "missing"),
          function(object, i) .zeroIsNA(object))

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "integer"),
          function(object, i) {
              el <- experiments(object)
              for (ii in i)
                  el[[ii]] <- zeroIsNA(el[[ii]])
              replaceAssay(object, el)
          })

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "numeric"),
          function(object, i) zeroIsNA(object, as.integer(i)))

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "character"),
          function(object, i) {
              if (any(! i %in% names(object)))
                  stop("subscript contains invalid names")
              zeroIsNA(object, which(names(object) %in% i))
          })


####---- infIsNA ----####


##' @exportMethod infIsNA
##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("SummarizedExperiment", "missing"),
          function(object, i) .infIsNA(object))

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "integer"),
          function(object, i) {
              el <- experiments(object)
              for (ii in i)
                  el[[ii]] <- infIsNA(el[[ii]])
              replaceAssay(object, el)
          })

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "numeric"),
          function(object, i) infIsNA(object, as.integer(i)))

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "character"),
          function(object, i) {
              if (any(! i %in% names(object)))
                  stop("subscript contains invalid names")
              infIsNA(object, which(names(object) %in% i))
          })


####---- nNA ----####


##' @exportMethod nNA
##' @rdname QFeatures-missing-data
setMethod("nNA", c("SummarizedExperiment", "missing"),
          function(object, i, addToObject = FALSE) .nNA(object, addToObject))

##' @rdname QFeatures-missing-data
setMethod("nNA", c("QFeatures", "integer"),
          function(object, i, addToObject = FALSE) .nNAi(object, i, addToObject))

##' @rdname QFeatures-missing-data
setMethod("nNA", c("QFeatures", "numeric"),
          function(object, i, addToObject = FALSE) .nNAi(object, as.integer(i), addToObject))

##' @rdname QFeatures-missing-data
setMethod("nNA", c("QFeatures", "character"),
          function(object, i, addToObject = FALSE) .nNAi(object, i, addToObject))


####---- filterNA ----####


##' @exportMethod filterNA
##' @rdname QFeatures-missing-data
setMethod("filterNA", "SummarizedExperiment",
          function(object, pNA = 0) {
              k <- .row_for_filterNA(assay(object), pNA)
              object[k, ]
          })

##' @rdname QFeatures-missing-data
setMethod("filterNA", "QFeatures",
          function(object, pNA = 0, i) {
              if (!is.character(i)) i <- names(object)[i]
              sel <- lapply(i, function(ii) {
                  .row_for_filterNA(assay(object[[ii]]), pNA)
              })
              names(sel) <- i
              object[sel, ]
          })
