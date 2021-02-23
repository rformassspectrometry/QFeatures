## ------------------
## Private functions
## ------------------

main_assay <- function(object)
    which.max(vapply(experiments(object),
                     nrow,
                     numeric(1)))

number_assays_in_se <- function(object)
    lengths(lapply(experiments(object), assays))


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
    al_names <- unname(vapply(object@assayLinks, "slot", "name",
                              FUN.VALUE = character(1)))
    ## An AssayLinks object is valid if the names of the node for all assays
    ## are contained in the assay names of the QFeatures object. The order
    ## matters!
    if (!all(al_names == n_exp))
        stop("@names not valid")
    al_from <- unname(unlist(lapply(object@assayLinks, "slot", "from")))
    ## An AssayLinks object is valid if the names of the parent assays for all
    ## assays are either NA (= root node) or contained in the assay names of the
    ## QFeatures object
    if (!all(is.na(al_from) | al_from %in% n_exp))
        stop("@from not valid")
    ## An AssayLinks object is valid if `hits` is not empty, unless 
    ## the assay is empty or it is the parent node
    hitsIsValid <- mapply(function(x, y){
        length(y@hits) != 0 || is.na(y@from) || nrow(x) == 0
    }, x = experiments(object), y = object@assayLinks)
    if (any(!hitsIsValid))
        stop("@hits are not valid")
    NULL
}

.unique_row_names <- function(object) {
    dup_row_names <- vapply(experiments(object),
                            function(x) anyDuplicated(rownames(x)),
                            numeric(1))
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
