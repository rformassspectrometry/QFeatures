##' @title Filter features based on their rowData
##'
##' @description
##'
##' The `filterFeatures` methods enables users to filter features
##' based on a variable in their `rowData`. The features matching the
##' filter will be returned as a new object of class `QFeatures`. The
##' filters can be provided as instances of class `AnnotationFilter`
##' (see below) or as formulas.
##'
##' @section The filtering procedure:
##'
##' `filterFeatures()` will go through each assay of the `QFeatures`
##' object and apply the filtering on the corresponding `rowData`.
##' Features that do not pass the filter condition are removed from
##' the assay. In some cases, one may want to filter for a variable
##' present in some assay, but not in other. There are two options:
##' either provide `keep = FALSE` to remove all features for those
##' assays (and thus leaving an empty assay), or provide `keep = TRUE`
##' to ignore filtering for those assays.
##'
##' Because features in a `QFeatures` object are linked between different
##' assays with `AssayLinks`, the links are automatically updated.
##' However, note that the function doesn't propagate the filter to parent
##' assays. For example, suppose a peptide assay with 4 peptides is
##' linked to a protein assay with 2 proteins (2 peptides mapped per
##' protein) and you apply `filterFeatures()`. All features pass the
##' filter except for one protein. The peptides mapped to that protein
##' will remain in the `QFeatures` object. If propagation of the
##' filtering rules to parent assay is desired, you may want to use
##' `x[i, , ]` instead (see the *Subsetting* section in `?QFeature`).
##'
##' @section Variable filters:
##'
##' The variable filters are filters as defined in the
##' [AnnotationFilter] package. In addition to the pre-defined filter,
##' users can arbitrarily set a field on which to operate. These
##' arbitrary filters operate either on a character variables (as
##' `CharacterVariableFilter` objects) or numerics (as
##' `NumericVariableFilters` objects), which can be created with the
##' `VariableFilter` constructor.
##'
##' @seealso The [QFeatures] man page for subsetting and the `QFeatures`
##'     vignette provides an extended example.
##'
##' @return An filtered `QFeature` object.
##'
##' @author Laurent Gatto
##'
##' @name QFeatures-filtering
##'
##' @rdname QFeatures-filtering
##'
##' @aliases filterFeatures filterFeatures,QFeatures,formula-method
##' @aliases filterFeatures,QFeatures,AnnotationFilter-method
##' @aliases CharacterVariableFilter NumericVariableFilter VariableFilter
##'
##' @examples
##'
##' ## ----------------------------------------
##' ## Creating character and numberic
##' ## variable filters
##' ## ----------------------------------------
##'
##' VariableFilter(field = "my_var",
##'                value = "value_to_keep",
##'                condition = "==")
##'
##' VariableFilter(field = "my_num_var",
##'                value = 0.05,
##'                condition = "<=")
##'
##' example(aggregateFeatures)
##'
##' ## ----------------------------------------------------------------
##' ## Filter all features that are associated to the Mitochondrion in
##' ## the location feature variable. This variable is present in all
##' ## assays.
##' ## ----------------------------------------------------------------
##'
##' ## using the forumla interface, exact mathc
##' filterFeatures(feat1, ~  location == "Mitochondrion")
##'
##' ## using the forumula intefrace, martial match
##' filterFeatures(feat1, ~startsWith(location, "Mito"))
##'
##' ## using a user-defined character filter
##' filterFeatures(feat1, VariableFilter("location", "Mitochondrion"))
##'
##' ## using a user-defined character filter with partial match
##' filterFeatures(feat1, VariableFilter("location", "Mito", "startsWith"))
##' filterFeatures(feat1, VariableFilter("location", "itochon", "contains"))
##'
##' ## ----------------------------------------------------------------
##' ## Filter all features that aren't marked as unknown (sub-cellular
##' ## location) in the feature variable
##' ## ----------------------------------------------------------------
##'
##' ## using a user-defined character filter
##' filterFeatures(feat1, VariableFilter("location", "unknown", condition = "!="))
##'
##' ## using the forumula interface
##' filterFeatures(feat1, ~ location != "unknown")
##'
##' ## ----------------------------------------------------------------
##' ## Filter features that have a p-values lower or equal to 0.03
##' ## ----------------------------------------------------------------
##'
##' ## using a user-defined numeric filter
##' filterFeatures(feat1, VariableFilter("pval", 0.03, "<="))
##'
##' ## using the formula interface
##' filterFeatures(feat1, ~ pval <= 0.03)
##'
##' ## you can also remove all p-values that are NA (if any)
##' filterFeatures(feat1, ~ !is.na(pval))
##'
##' ## ----------------------------------------------------------------
##' ## Negative control - filtering for an non-existing markers value,
##' ## returning empty results.
##' ## ----------------------------------------------------------------
##'
##' filterFeatures(feat1, VariableFilter("location", "not"))
##'
##' filterFeatures(feat1, ~ location == "not")
##'
##' ## ----------------------------------------------------------------
##' ## Filtering for a  missing feature variable. The outcome is controled
##' ## by keep
##' ## ----------------------------------------------------------------
##' data(feat2)
##'
##' filterFeatures(feat2, ~ y < 0)
##'
##' filterFeatures(feat2, ~ y < 0, keep = TRUE)
##'
##' ## ----------------------------------------------------------------
##' ## Example with missing values
##' ## ----------------------------------------------------------------
##'
##' data(feat1)
##' rowData(feat1[[1]])[1, "location"] <- NA
##' rowData(feat1[[1]])
##'
##' ## The row with the NA is not removed
##' rowData(filterFeatures(feat1, ~ location == "Mitochondrion")[[1]])
##' rowData(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = FALSE)[[1]])
##'
##' ## The row with the NA is removed
##' rowData(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = TRUE)[[1]])
##'
##' ## Note that in situations with missing values, it is possible to
##' ## use the `%in%` operator or filter missing values out
##' ## explicitly.
##'
##' rowData(filterFeatures(feat1, ~ location %in% "Mitochondrion")[[1]])
##' rowData(filterFeatures(feat1, ~ location %in% c(NA, "Mitochondrion"))[[1]])
##'
##' ## Explicit handling
##' filterFeatures(feat1, ~ !is.na(location) & location == "Mitochondrion")
##'
##' ## Using the pipe operator
##' feat1 |>
##'    filterFeatures( ~ !is.na(location)) |>
##'    filterFeatures( ~ location == "Mitochondrion")
NULL



##' @import AnnotationFilter
##' @exportClass CharacterVariableFilter
##' @rdname QFeatures-filtering
setClass("CharacterVariableFilter", contains = "CharacterFilter")

##' @exportClass NumericVariableFilter
##' @rdname QFeatures-filtering
setClass("NumericVariableFilter", contains = "DoubleFilter")


##' @param field `character(1)` refering to the name of the variable
##'     to apply the filter on.
##'
##' @param value `character()` or `integer()` value for the
##'     `CharacterVariableFilter` and `NumericVariableFilter` filters
##'     respectively.
##'
##' @param condition `character(1)` defining the condition to be used in
##'     the filter. For `NumericVariableFilter`, one of `"=="`,
##'     `"!="`, `">"`, `"<"`, `">="` or `"<="`. For
##'     `CharacterVariableFilter`, one of `"=="`, `"!="`,
##'     `"startsWith"`, `"endsWith"` or `"contains"`. Default
##'     condition is `"=="`.
##'
##' @param not `logical(1)` indicating whether the filtering should be negated
##'     or not. `TRUE` indicates is negated (!). `FALSE` indicates not negated.
##'     Default `not` is `FALSE`, so no negation.
##'
##' @export VariableFilter
##' @rdname QFeatures-filtering
VariableFilter <- function(field,
                           value,
                           condition = "==",
                           not = FALSE) {
    if (is.numeric(value))
        new("NumericVariableFilter",
            field = as.character(field),
            value = value,
            condition = condition,
            not = not)
    else if (is.character(value))
        new("CharacterVariableFilter",
            field = as.character(field),
            value = value,
            condition = condition,
            not = not)
    else
        stop("Value type undefined.")
}



##' @param object An instance of class [QFeatures].
##'
##' @param filter Either an instance of class [AnnotationFilter] or a
##'     formula.
##'
##' @param i A numeric, logical or character vector pointing to the
##'     assay(s) to be filtered.
##'
##' @param na.rm `logical(1)` indicating whether missing values should
##'     be removed. Default is `FALSE`.
##'
##' @param keep `logical(1)` indicating whether to keep the features
##'     of assays for which at least one of the filtering variables are
##'     missing in the rowData. When `FALSE` (default), all such assay
##'     will contain 0 features; when `TRUE`, the assays are untouched.
##'
##' @param ... Additional parameters. Currently ignored.
##'
##' @exportMethod filterFeatures
##'
##' @importMethodsFrom ProtGenerics filterFeatures
##'
##' @rdname QFeatures-filtering
setMethod("filterFeatures",
          c("QFeatures", "AnnotationFilter"),
          function(object, filter, i, na.rm = FALSE,
                   keep = FALSE, ...)
              filterFeaturesWithAnnotationFilter(object, filter, i,
                                                 na.rm, keep, ...))
##' @rdname QFeatures-filtering
setMethod("filterFeatures",
          c("QFeatures", "formula"),
          function(object, filter, i, na.rm = FALSE,
                   keep = FALSE, ...)
              filterFeaturesWithFormula(object, filter, i,
                                        na.rm, keep, ...)
          )

##' @importFrom BiocGenerics do.call
filterFeaturesWithAnnotationFilter <- function(object, filter, i,
                                               na.rm, keep, ...) {
    ## Check the index
    i <- if (missing(i)) names(object) else .normIndex(object, i)

    ## Check the filtering variables
    vars <- field(filter)
    isPresent <- .checkFilterVariables(rowData(object), vars)

    ## Apply the filter
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      if (field(filter) %in% names(x))
                          do.call(condition(filter),
                                  list(x[, field(filter)],
                                       value(filter)))
                      else
                          rep(FALSE, nrow(x))
                  })
    ## Take missing data into account
    sel <- lapply(sel, function(x) {
        x[is.na(x)] <- !na.rm
        x
    })

    ## Apply the VariableFilter inversion
    if (not(filter)) sel <- lapply(sel, "!")

    ## If required, keep lost assays
    if (keep) sel <- .keepLostAssays(sel, isPresent)

    ## Reset the filter for assays not selected by i
    for (ii in names(sel)[!names(sel) %in% i]) {
        sel[[ii]][] <- TRUE
    }

    object[sel, , ]
}


##' @importFrom lazyeval f_eval
filterFeaturesWithFormula <- function(object, filter, i,
                                      na.rm, keep, ...) {
    ## Check the index
    i <- if (missing(i)) names(object) else .normIndex(object, i)

    ## Check the filtering variables
    vars <- all.vars(filter)
    isPresent <- .checkFilterVariables(rowData(object), vars)

    ## Apply the filter
    sel <- lapply(experiments(object),
                  function(exp) {
                      x <- rowData(exp)
                      tryCatch(lazyeval::f_eval(filter, data = as.list(x)),
                               error = function(e) rep(FALSE, nrow(x)))
                  })

    ## Take missing data into account
    sel <- lapply(sel, function(x) {
        x[is.na(x)] <- !na.rm
        x
    })

    ## If required, keep lost assays
    if (keep) sel <- .keepLostAssays(sel, isPresent)

    ## Reset the filter for assays not selected by i
    for (ii in names(sel)[!names(sel) %in% i]) {
        sel[[ii]][] <- TRUE
    }

    object[sel, , ]
}

##' Internal function that checks whether the variables queried by the
##' filter are available in a given a List of DataFrame objects, each
##' element representing the rowData associated with an assay in a
##' QFeatures object.
##'
##' If one or more variables are missing from *all* assays, an error
##' is thrown because no filtering can be applied. The function will
##' ignore variables that are stored in parent environments as these
##' are used as filtering values rather than variables. The function
##' prints the search result to the console as well, to inform the
##' user of how many assays contain the desired filtering variable(s),
##' and which assays are ignored due to missing filtering
##' variables(s).
##'
##' @param rowdata `List` of `DataFrame` objects reprensting the
##'     rowData of a `QFeatures` object. Each element is related to an
##'     assay.
##'
##' @param vars `character()` with one or more filter variable names.
##'     The function checks their presence. Variables present in
##'     parent environment are ignored.
##'
##' @param parent.frame `logical(1)` should variables that are present
##'     in `parent.frame(4)` be removed. Default is `TRUE`. Set to
##'     `FALSE` when using `VariableFilter()`, as only filter values
##'     are passed.
##'
##' @return A `matrix` where rows represent filter variables
##'     (excluding variables found in the parent environment) and
##'     columns represent assays. Each element is a logical indicating
##'     whether a given variable was found in the rowData associated
##'     to a given assay (`TRUE`) or not (`FALSE`).
##'
##' @noRd
.checkFilterVariables <- function(rowdata, vars) {
    ## Ignore variables from the user environment. We search for
    ## variables to omit from the check in the 4th parent environment
    ## (may not always be .GlobalEnv). Here is a "traceback" counter:
    ## 0 in .checkFilterVariables()
    ## 1 in FilterFeaturesWithFormula()
    ## 2 in .local()
    ## 3 in filterFeatures()
    ## 4 in environment the function was called
    ##
    ## This is needed for when the value is a variable itself, such as
    ## below, because we don't want to search for target in the
    ## rowData.
    ##
    ## target <- "location"
    ## filterFeatures(feat1, ~  location == target)
    ##
    ## BUT this breaks if location exists in the working env
    ## (described also in issue #208)
    ##
    ## location <- 1
    ## filterFeatures(feat1, ~  location == "Mitchondrion")
    ## filterFeatures(feat1, ~  location == target)
    ##
    ## The number of variables (that we want to keep, vs their values
    ## (that we don't want) isn't necessarily 1, as shown in:
    ## filterFeatures(feat1, ~ pval <= 0.03 & grepl("Mito", location))
    ##
    v1 <- vars[1]
    vars <- vars[!vars %in% ls(envir = parent.frame(4))]
    ## vars[1], should always be a proper var (unless there's a
    ## typo). Add it back to avoid removing all vars.
    vars <- unique(c(vars, v1))
    if (!length(vars)) ## this should never happen anymore
        stop("No filter variables left.")
    ## get in which assays each variable comes from
    out <- sapply(colnames(rowdata), function(rdn) vars %in% rdn)
    if (!is.array(out)) out <- t(out)
    rownames(out) <- vars
    ## Throw an error if variables are missing from all assays
    err <- rowSums(out) == 0
    if (any(err)) stop("'", paste(names(err)[err], collapse = "', '"),
                       "' is/are absent from all rowData.")
    ## Get the assays for which one of the variables is absent in the rd
    absent <- colnames(out)[colSums(!out) > 0]
    ## Print a message for each variable
    msg <- sapply(vars, function(var) {
        x <- sum(out[var, ])
        paste0("'", var, "' found in ", x, " out of ", length(rowdata),
               " assay(s)\n")
    })
    if (length(absent) > 0)
        msg <- c(msg, "No filter applied to the following assay(s) because ",
                 "one or more filtering variables are missing ",
                 "in the rowData: ", paste0(absent, collapse = ", "), ".\n",
                 "You can control whether to remove or keep the features ",
                 "using the 'keep' argument (see '?filterFeature').")
    message(msg)
    ## Return the presence matrix
    out
}

## Internal function to keep the all rows for assays that contain no
## filter variable among their rowData.
## @param x A list of logical(). Each element of the list represents
##     an assay and contains a vector of length equal to the number of
##     rows (features) of the associated assay. This vector contains
##     logicals indicating whether the corresponding row (feature) is
##     kept (TRUE) or removed (FALSE).
## @param mat A matrix as return by .checkFilterVariables().
## @return A modified version of x, where elements that contain a
##     vector with all FALSE were changed to all TRUE.
.keepLostAssays <- function(x, mat) {
    ## Get assays that have at least one missing filter variable
    ## Keep all features for those assays
    absent <- colnames(mat)[colSums(!mat) > 0]
    for (i in absent) x[[i]][] <- TRUE
    x
}

## Internal function called by `filterFeaturesWithAnnotationFilter` when
## `condition` is `"contains"`
contains <- function(x, value) {
    ## Replace regex special character by regular character matching
    value <- gsub('([[:punct:]])', '\\[\\1\\]', value)
    ## Return whether elements in x contain value or not
    grepl(pattern = value, x = x)
}
