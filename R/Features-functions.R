## ------------------
## Private functions
## ------------------

main_assay <- function(object)
    which.max(sapply(object, nrow))

## ----------------------------
## Internal validity functions
## ----------------------------

##' @importFrom methods slot
.valid_Features_indices <- function(object) {
    ## ids <- get_featureSet_ids(object)
    ## if (anyNA(ids)) 
    ##     stop("Indices musn't be NA")
    ## if (anyDuplicated(ids))
    ##     stop("Duplicated indices")
    NULL
}

.valid_Features <- function(object) {
    .valid_Features_indices(object)
}

setValidity("Features", .valid_Features)
