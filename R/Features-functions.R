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

.valid_assay_links <- function(object) {
    n_exp <- names(experiments(object))
    al_names <- unname(sapply(object@assayLinks, "slot", "name"))
    ## An AssayLinks object is valid if the names of the node for all assays 
    ## are contained in the assay names of the Features object
    if (!all(al_names %in% n_exp))
        stop("@names not valid")
    al_from <- unname(unlist(sapply(object@assayLinks, "slot", "from")))
    ## An AssayLinks object is valid if the names of the parent assays for all 
    ## assays are either NA (= root node) or contained in the assay names of the 
    ## Features object
    if (!all(is.na(al_from) | al_from %in% n_exp))
        stop("@from not valid")
    NULL
}

.valid_assays <- function(object) {
    if (any(!sapply(experiments(object), inherits, "SummarizedExperiment")))
        stop("Invalid assay. All assays must inherit from 'SummarizedExperiment' class.")
}

.valid_Features <- function(object) {
    .valid_Features_indices(object)
    .valid_assay_links(object)
    .valid_assays(object)
}

setValidity("Features", .valid_Features)
