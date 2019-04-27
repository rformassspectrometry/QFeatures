## ----------------------
## id and from management
## ----------------------

##' Returns the [FeatureSet]s `id` slots.
##'
##' @param object A Features object.
##' @return The [FeatureSet]s `id` slots.
##' @md
##' @rdname FeatureSet-links
get_featureSet_ids <- function(object) 
    sapply(object, slot, "id")


##' This function sets the [FeatureSet]'s `id` slot to
##' `1:length(object)`. The [FeatureSet] with the larges number of
##' features gets `id` 1, next one 2, ... and the smallest one get
##' `length(object)`.
##' 
##' @return A Features object with FeatureSet element ids set to
##'     `1:length(object)`.
##' @md
##' @rdname FeatureSet-links
set_featureSet_ids <- function(object) {
    if (isEmpty(object))
        return(object)
    ln <- length(object)
    new_ids <- seq_len(ln)
    o <- order(dims(object)[1, ], decreasing = TRUE)
    new_ids <- new_ids[o]
    for (i in seq_len(ln))
        object@listData[[i]]@id <- new_ids[i]
    if (validObject(object))
        object
}

##' Returns the next identifier.
##'
##' @return The next `id` in the [Features], as an `integer(1)`.
get_next_featureSet_id <- function(object) 
    as.integer(max(get_featureSet_ids(object)) + 1)
