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
    n_alk <- names(object@assayLinks)
    if (any(n_exp != n_alk))
        stop("Assay and link names don't match")
    al_names <- unname(sapply(object@assayLinks, "slot", "name"))
    if (!is.na(al_names) && !all(al_names %in% n_exp))
        stop("@names not valid")
    al_from <- unname(sapply(object@assayLinks, "slot", "from"))
    if (!is.na(al_from) && !all(al_from %in% n_exp))
        stop("@from not valid")
    NULL
}

.valid_Features <- function(object) {
    .valid_Features_indices(object)
    .valid_assay_links(object)
}

setValidity("Features", .valid_Features)
