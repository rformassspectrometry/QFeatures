## ------------------
## Private functions
## ------------------

main_assay <- function(object)
    which.max(sapply(experiments(object), nrow))

number_assays_in_se <- function(object) {
    lengths(sapply(experiments(object), assays))    
}

## ----------------------------
## Internal validity functions
## ----------------------------

##' @importFrom methods slot
.valid_Features_indices <- function(object) {
    if (!isEmpty(object) && !identical(names(object), names(object@assayLinks)))
        stop("Assay links names are wrong.")    
    NULL
}

.valid_Features <- function(object) {
    .valid_Features_indices(object)
}

setValidity("Features", .valid_Features)
