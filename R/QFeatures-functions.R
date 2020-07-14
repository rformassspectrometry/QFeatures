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
.valid_QFeatures_indices <- function(object) {
    if (!isEmpty(object) && !identical(names(object), names(object@assayLinks)))
        stop("Assay links names are wrong.")
    NULL
}

.valid_assay_links <- function(object) {
    n_exp <- names(experiments(object))
    al_names <- unname(sapply(object@assayLinks, "slot", "name"))
    ## An AssayLinks object is valid if the names of the node for all assays
    ## are contained in the assay names of the QFeatures object
    if (!all(al_names %in% n_exp))
        stop("@names not valid")
    al_from <- unname(unlist(sapply(object@assayLinks, "slot", "from")))
    ## An AssayLinks object is valid if the names of the parent assays for all
    ## assays are either NA (= root node) or contained in the assay names of the
    ## QFeatures object
    if (!all(is.na(al_from) | al_from %in% n_exp))
        stop("@from not valid")
    NULL
}

.unique_row_names <- function(object) {
    dup_row_names <- sapply(experiments(object),
                            function(x) anyDuplicated(rownames(x)))
    if (any(dup_row_names != 0))
        stop("Assay(s) ", paste(which(dup_row_names != 0), collapse = ", "),
             " has/have duplicated row names.")
    NULL
}

.valid_QFeatures <- function(object) {
    .valid_QFeatures_indices(object)
    .valid_assay_links(object)
    .unique_row_names(object)
}

setValidity("QFeatures", .valid_QFeatures)
