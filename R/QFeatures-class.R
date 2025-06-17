##' @title Quantitative MS QFeatures
##'
##' @description
##'
##' Conceptually, a `QFeatures` object holds a set of *assays*, each
##' composed of a `matrix` (or `array`) containing quantitative data
##' and row annotations (meta-data).  The number and the names of the
##' columns (samples) must always be the same across the assays, but
##' the number and the names of the rows (features) can vary. The
##' assays are typically defined as `SummarizedExperiment` objects. In
##' addition, a `QFeatures` object also uses a single `DataFrame` to
##' annotate the samples (columns) represented in all the matrices.
##'
##' The `QFeatures` class extends the
##' [MultiAssayExperiment::MultiAssayExperiment] and inherits all
##' the functionality of the
##' [MultiAssayExperiment::MultiAssayExperiment] class.
##'
##' A typical use case for such `QFeatures` object is to represent
##' quantitative proteomics (or metabolomics) data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values. The largest assay (the one with the
##' highest number of features, PSMs in the example above) is
##' considered the main assay.
##'
##' The recommended way to create `QFeatures` objects is the use the
##' [readQFeatures()] function, that creates an instance from tabular
##' data. The `QFeatures` constructor can be used to create objects
##' from their bare parts.  It is the user's responsability to make
##' sure that these match the class validity requirements.
##'
##' @section Constructors:
##'
##' - `QFeatures(..., assayLinks)` allows the manual construction of
##'   objects. It is the user's responsability to make sure these
##'   comply. The arguments in `...` are those documented in
##'   [MultiAssayExperiment::MultiAssayExperiment()]. For details
##'   about `assayLinks`, see [AssayLinks]. An example is shown below.
##'
##' - The [readQFeatures()] function constructs a `QFeatures` object
##'   from text-based spreadsheet or a `data.frame` used to generate
##'   an assay. See the function manual page for details and an
##'   example.
##'
##' @section Accessors:
##'
##' - The `QFeatures` class extends the
##'   [MultiAssayExperiment::MultiAssayExperiment] class and inherits
##'   all its accessors and replacement methods.
##'
##' - The `rowData` method returns a `DataFrameList` containing the
##'   `rowData` for each assay of the `QFeatures` object. On the other
##'   hand, `rowData` can be modified using `rowData(x) <- value`,
##'   where `value` is a list of tables that can be coerced to `DFrame`
##'   tables. The names of `value` point to the assays for
##'   which the `rowData` must be replaced. The column names of each
##'   table are used to replace the data in the existing `rowData`. If
##'   the column name does not exist, a new column is added to the
##'   `rowData`.
##'
##' - The `rbindRowData` functions returns a `DFrame` table that
##'   contains the row binded `rowData` tables from the selected
##'   assays. In this context, `i` is a `character()`, `integer()` or
##'   `logical()` object for subsetting assays. Only rowData variables
##'   that are common to all assays are kept.
##'
##' - The `rowDataNames` accessor returns a list with the `rowData`
##'   variable names.
##'
##' - The [longForm()] accessor takes a `QFeatures` instance and returns it in a
##'   long *tidy* `DataFrame`, where each quantitative value is reported on a
##'   separate line.
##'
##' @section Adding, removing and replacing assays:
##'
##' - The [aggregateFeatures()] function creates a new assay by
##'   aggregating features of an existing assay.
##'
##' - `addAssay(x, y, name, assayLinks)`: Adds one or more
##'   new assay(s) `y` to the `QFeatures` instance `x`. `name`
##'   is a `character(1)` naming the assay if only one assay is
##'   provided, and is ignored if `y` is a list of assays. `assayLinks`
##'   is an optional [AssayLinks]. The `colData(y)` is
##'   automatically added to `colData(x)` by matching sample
##'   names, that is `colnames(y)`. If the samples are not present in
##'   `x`, the rows of `colData(x)` are extended to account for the
##'   new samples. Be aware that conflicting information between the
##'   `colData(y)` and the `colData(x)` will result in an
##'   error.
##'
##' - `removeAssay(x, i)`: Removes one or more assay(s) from the
##'   `QFeatures` instance `x`. In this context, `i` is a `character()`,
##'   `integer()` or `logical()` that indicates which assay(s) to
##'   remove.
##'
##' - `replaceAssay(x, y, i)`: Replaces one or more
##'   assay(s) from the `QFeatures` instance `x`. In this context, `i`
##'   is a `character()`, `integer()` or `logical()` that indicates
##'   which assay(s) to replace. The `AssayLinks` from or to
##'   any replaced assays are automatically removed, unless the
##'   replacement has the same dimension names (columns and row, order
##'   agnostic). Be aware that conflicting information between
##'   `colData(y)` and `colData(x)` will result in an error.
##'
##' - `x[[i]] <- value`: a generic method for adding (when `i` is not
##'   in `names(x)`), removing (when `value` is null) or replacing (when
##'   `i` is in `names(x)`). Note that the arguments `j` and `...` from
##'   the S4 replacement method signature are not allowed.
##'
##' @section Subsetting:
##'
##' - QFeatures object can be subset using the `x[i, j, k, drop = TRUE]`
##'   paradigm. In this context, `i` is a `character()`, `integer()`,
##'   `logical()` or `GRanges()` object for subsetting by rows. See
##'   the argument descriptions for details on the remaining arguments.
##'
##' - The [subsetByFeature()] function can be used to subset a
##'   `QFeatures` object using one or multiple feature names that will
##'   be matched across different assays, taking the aggregation
##'   relation between assays.
##'
##' - The `selectRowData(x, rowvars)` function can be used to
##'   select a limited number of `rowData` columns of interest named
##'   in `rowvars` in the `x` instance of class `QFeatures`. All other
##'   variables than `rowvars` will be dropped. In case an element in
##'   `rowvars` isn't found in any `rowData` variable, a message is
##'   printed.
##'
##' - The `dropEmptyAssays(object, dims)` function removes empty
##'   assays from a `QFeatures`. Empty assays are defined as having 0
##'   rows and/or 0 columns, as defined by the `dims` argument.
##'
##' @param object An instance of class [QFeatures].
##'
##' @param x An instance of class [QFeatures].
##'
##' @param i An indexing vector. See the corresponding section in the
##'     documentation for more details.
##'
##' @param j `character()`, `logical()`, or `numeric()` vector for
##'     subsetting by `colData` rows.
##'
##' @param k `character()`, `logical()`, or `numeric()` vector for
##'     subsetting by assays
##'
##' @param value The values to use as a replacement. See the
##'     corresponding section in the documentation for more details.
##'
##' @param drop logical (default `TRUE`) whether to drop empty assay
##'     elements in the `ExperimentList`.
##'
##' @param rowvars A `character()` with the names of the `rowData`
##'     variables (columns) to retain in any assay.
##'
##' @param ... See `MultiAssayExperiment` for details. For `plot`,
##'     further arguments passed to `igraph::plot.igraph`.
##'
##'
##'
##' @seealso
##'
##' - The [readQFeatures()] constructor and the [aggregateFeatures()]
##'   function. The *QFeatures* vignette provides an extended example.
##'
##' - The [QFeatures-filtering] manual page demonstrates how to filter
##'   features based on their rowData.
##'
##' - The [missing-data] manual page to manage missing values in
##'   `QFeatures` objects.
##'
##' - The [QFeatures-processing] and [aggregateFeatures()] manual pages
##'   and *Processing* vignette describe common quantitative data
##'   processing methods using in quantitative proteomics.
##'
##' @import MultiAssayExperiment ProtGenerics
##'
##' @name QFeatures
##'
##' @aliases QFeatures QFeatures-class class:QFeatures
##' @aliases addAssay
##' @aliases dims,QFeatures-method show,QFeatures-method
##' @aliases [,QFeatures,ANY,ANY,ANY-method [,QFeatures,character,ANY,ANY-method
##' @aliases rowDataNames selectRowData
##'
##' @rdname QFeatures-class
##'
##' @exportClass QFeatures
##'
##' @author Laurent Gatto
##'
##' @examples
##' ## ------------------------
##' ## An empty QFeatures object
##' ## ------------------------
##'
##' QFeatures()
##'
##' ## -----------------------------------
##' ## Creating a QFeatures object manually
##' ## -----------------------------------
##'
##' ## two assays (matrices) with matching column names
##' m1 <- matrix(1:40, ncol = 4)
##' m2 <- matrix(1:16, ncol = 4)
##' sample_names <- paste0("S", 1:4)
##' colnames(m1) <- colnames(m2) <- sample_names
##' rownames(m1) <- letters[1:10]
##' rownames(m2) <- letters[1:4]
##'
##' ## two corresponding feature metadata with appropriate row names
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = rownames(m1))
##' df2 <- DataFrame(row.names = rownames(m2))
##'
##' (se1 <- SummarizedExperiment(m1, df1))
##' (se2 <- SummarizedExperiment(m2, df2))
##'
##' ## Sample annotation (colData)
##' cd <- DataFrame(Var1 = rnorm(4),
##'                 Var2 = LETTERS[1:4],
##'                 row.names = sample_names)
##'
##' el <- list(assay1 = se1, assay2 = se2)
##' fts1 <- QFeatures(el, colData = cd)
##' fts1
##' fts1[[1]]
##' fts1[["assay1"]]
##'
##' ## Rename assay
##' names(fts1) <- c("se1", "se2")
##'
##' ## Add an assay
##' fts1 <- addAssay(fts1, se1[1:2, ], name = "se3")
##'
##' ## Get the assays feature metadata
##' rowData(fts1)
##'
##' ## Keep only the Fa variable
##' selectRowData(fts1, rowvars = "Fa")
##'
##' ## -----------------------------------
##' ## See ?readQFeatures to create a
##' ## QFeatures object from a data.frame
##' ## or spreadsheet.
##' ## -----------------------------------
##'
NULL


## ----------------------------------
## QFeatures Class ChangeLog
##
## Version 0.1:
##  - Contains MatchedAssayExperiment
## Version 0.2:
##  - Contains MultiAssayExperiment (see issue 46)
## Version 0.3:
##  - Rename to QFeatures (see issue 89)

setClass("QFeatures",
    contains = "MultiAssayExperiment",
    slots = c(
        version = "character",
        assayLinks = "AssayLinks"
    ),
    prototype = prototype(
        version = "0.3"
    )
)

##' @export
##'
##' @rdname QFeatures-class
##'
##' @param assayLinks An optional [AssayLinks] object.
QFeatures <- function(..., assayLinks = NULL) {
    ans <- MultiAssayExperiment(...)
    if (isEmpty(ans)) {
        assayLinks <- AssayLinks()
    } else {
        if (is.null(assayLinks)) {
            assayLinks <- AssayLinks(names = names(ans))
        }
    }
    new("QFeatures",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = assayLinks
    )
}


##' @rdname QFeatures-class
##'
##' @exportMethod show
setMethod(
    "show", "QFeatures",
    function(object) {
        # suppress messages in case of implicit QFeatures type
        type <- suppressMessages(getQFeaturesType(object))
        if (isEmpty(object)) {
            cat(sprintf(
                "An empty instance of class %s (type: %s)\n",
                class(object), type
            ))
            return(NULL)
        }
        n <- length(object)
        cat(sprintf(
            "An instance of class %s (type: %s) with %d set%s:\n",
            class(object), type, n, ifelse(n == 1, "", "s")
        ))
        el <- experiments(object)
        o_class <- class(el)
        elem_cl <- vapply(el, class, character(1L))
        o_len <- length(el)
        o_names <- names(el)
        featdim <- vapply(el, FUN = function(obj) {
            dim(obj)[1]
        }, FUN.VALUE = integer(1L))
        sampdim <- vapply(el, FUN = function(obj) {
            dim(obj)[2]
        }, FUN.VALUE = integer(1L))
        if (n <= 7) {
            cat(sprintf(
                "\n [%i] %s: %s with %s rows and %s columns",
                seq(o_len), o_names, elem_cl,
                featdim, sampdim
            ), "\n")
        } else {
            cat(sprintf(
                "\n [%i] %s: %s with %s rows and %s columns",
                seq(o_len)[1:3], o_names[1:3], elem_cl[1:3],
                featdim[1:3], sampdim[1:3]
            ), "\n")
            cat(" ...")
            cat(sprintf(
                "\n [%i] %s: %s with %s rows and %s columns",
                seq(o_len)[(n - 2):n], o_names[(n - 2):n], elem_cl[(n - 2):n],
                featdim[(n - 2):n], sampdim[(n - 2):n]
            ), "\n")
        }
    }
)


## Function that creates a plotly network graph from an igraph object
## @param graph An i graph object.
## @param coords A n vertices by 2 matrix with the coordinates of the
##     nodes.
##' @importFrom igraph get.edgelist layout_as_tree plot.igraph add_edges
##' @importFrom grDevices rgb
##'
.plotlyGraph <- function(graph, coords) {
    stopifnot(inherits(graph, "igraph"))
    ## Initialize plotly
    pl <- plotly::plot_ly()
    ## Add edges
    el <- get.edgelist(graph)
    if (nrow(el) > 0) {
        edge_coords <- sapply(1:nrow(el), function(i) {
            edge_coord <- c(
                coords[el[i, 1], 1],
                coords[el[i, 1], 2],
                coords[el[i, 2], 1],
                coords[el[i, 2], 2]
            )
        })
        pl <- plotly::add_segments(pl,
            x = edge_coords[1, ],
            y = edge_coords[2, ],
            xend = edge_coords[3, ],
            yend = edge_coords[4, ],
            line = list(
                color = "grey", width = 0.3,
                showarrow = TRUE
            )
        )
    }
    ## Add nodes
    pl <- plotly::add_markers(pl,
        x = coords[, 1], y = coords[, 2],
        marker = list(
            color = rgb(0.8, 0.8, 0.8),
            size = 40
        ),
        text = names(V(graph)),
        hoverinfo = "text"
    )
    ## Add labels
    pl <- plotly::add_text(pl,
        x = coords[, 1], y = coords[, 2],
        hoverinfo = "text",
        text = names(V(graph))
    )

    ## Edit plot plot
    axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
    plotly::layout(pl,
        xaxis = axis,
        yaxis = axis,
        showlegend = FALSE
    )
}

## Offset on the coordinates for better rendering when many assays
## have to be drawn
.offsetNodes <- function(coords) {
    lev <- coords[, 2]
    nlev <- length(unique(lev))
    ## Compute the step between levels, if only a single level the step=1
    step <- ifelse(nlev > 1, max(coords[, 2]) / (nlev - 1), 1)
    ## Custom margin between nodes so that the nodes are interleaved
    ## (pattern repeats every 6 node)
    mar <- c(1, 3, 5, 2, 4, 6) / 10
    for (i in unique(lev)) {
        sel <- lev == i
        ## Center x
        coords[sel, 1] <- coords[sel, 1] - mean(coords[sel, 1])
        ## Offset y
        mari <- rep(mar, length.out = sum(sel))
        mari <- mari - mean(mari)
        offset <- mari * step
        coords[sel, 2] <- coords[sel, 2] + offset
    }
    coords
}

##' @rdname QFeatures-class
##'
##' @param interactive A `logical(1)`. If `TRUE`, an interactive graph
##'     is generated using `plotly`. Else, a static plot using `igraph`
##'     is generated. We recommend interactive exploration when the
##'     `QFeatures` object contains more than 50 assays.
##'
##' @importFrom igraph make_graph layout_as_tree plot.igraph add_edges V
##' @export
plot.QFeatures <- function(x, interactive = FALSE, ...) {
    ## Check arguments
    if (!interactive & length(x) > 50) {
        warning(
            "The QFeatures object contains many assays. You may ",
            "want to consider creating an interactive plot (set ",
            "'interactive = TRUE')"
        )
    }
    ## Create the network graph
    graph <- make_graph(
        edges = character(0),
        isolates = names(x)
    )
    ## Add the edges = links between assays
    roots <- c()
    for (child in names(x)) {
        parents <- assayLink(x, child)@from
        for (parent in parents) {
            if (!is.na(parent)) {
                graph <- add_edges(graph, c(parent, child))
            } else {
                roots <- c(roots, child)
            }
        }
    }
    ## Tree layout
    coords <- layout_as_tree(graph, root = roots)
    coords <- .offsetNodes(coords)
    rownames(coords) <- names(V(graph))
    ## Perform plotting
    if (!interactive) {
        plot.igraph(graph, layout = coords, ...)
        return(invisible(NULL))
    } else {
        return(.plotlyGraph(graph, coords))
    }
}


## Function that prunes a `Hits` object from an `AssayLink` object,
## making sure that the `Hits` object is still valid with respect to
## a parent assay and its corresponding assay (self). The validity is
## ensured by removing missing features.
##
## @param hits A `Hits` object
## @param parent A `SummarizedExperiment` object or any object that
##     inherits from it. This is the assay that `hits` links from.
## @param self A `SummarizedExperiment` object or any object that
##     inherits from it. This is the assay that `hits` links to.
##
.pruneHits <- function(hits, parent, self) {
    ## Get the feature names in the parent and self assay
    featnParent <- rownames(parent)
    featnSelf <- rownames(self)
    ## Check which links are still in parent and self
    inParent <- mcols(hits)$names_from %in% featnParent
    inSelf <- mcols(hits)$names_to %in% featnSelf
    ## Remove lost feature links
    hits[inParent & inSelf, ]
}

## Function that prunes an `AssayLink` of a `QFeatures` object,
## making sure that the `AssayLink` object is still valid with respect
## to a given `QFeatures` object.
##
## @param al An `AssayLink` object
## @param object A `QFeatures` object. `al` will be adapted so that it
##     becomes valid when contained in `object`.
##
.pruneAssayLink <- function(al, object) {
    ## Identify lost assays that need to be pruned
    lost <- !al@from %in% names(object)
    ## Prune the links to assays in `@name`. When an assay is lost,
    ## the links to that assay are removed.
    if (all(lost)) { ## If all parent assays are lost, return an
        ## empty AssayLink
        al <- AssayLink(name = al@name)
        return(al)
    } else if (any(lost)) { ## If some parents are lost (only in case
        ## of multiple parents), remove the link to the lost parent(s)
        ## (`@from`) and subset the corresponding `Hits` object(s)
        al@from <- al@from[!lost]
        al@hits <- al@hits[!lost]
    }
    ## Prune the links to features in `@hits`. Even when an assay is
    ## not lost, some of its features might be lost and the
    ## corresponding `Hits` must be adapted.
    if (inherits(al@hits, "List")) { ## If the AssayLink contains a
        ## HitsList object, iterate through each `Hits` element.
        al@hits <- mendoapply(function(hits, parent) {
            .pruneHits(hits, parent, self = object[[al@name]])
        }, hits = al@hits, parent = experiments(object)[al@from])
    } else { ## If the AssayLink contains a single Hits object
        al@hits <- .pruneHits(
            hits = al@hits,
            parent = object[[al@from]],
            self = object[[al@name]]
        )
    }
    al
}

## Function that prunes the `AssayLinks` of a `QFeatures` object,
## making sure that the `AssayLinks` object is still valid after
## `QFeatures` subsetting
##
## @param object A `QFeatures` object
##
.pruneAssayLinks <- function(object) {
    ## Subset the AssayLinks
    object@assayLinks <- object@assayLinks[names(object)]
    ## Removed lost links in each AssayLink object
    object@assayLinks <- endoapply(object@assayLinks,
        .pruneAssayLink,
        object = object
    )
    ## Check new AssaLinks are valid
    .validAssayLinks(object)
    object
}

##' @rdname QFeatures-class
##'
##' @importFrom methods callNextMethod
##'
##' @exportMethod [
setMethod(
    "[", c("QFeatures", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        ## Subset the assays
        ans <- callNextMethod(x, i, j, ..., drop)

        ## Prune the AssayLinks so that the `QFeatures` object
        ## remains valid
        .pruneAssayLinks(ans)
    }
)

##' @rdname QFeatures-class
setMethod(
    "[", c("QFeatures", "character", "ANY", "ANY"),
    function(x, i, j, k, ..., drop = TRUE) {
        if (missing(j)) j <- TRUE
        if (missing(k)) k <- TRUE
        subsetByFeature(x, i)[, j, k]
    }
)


##' @rdname QFeatures-class
##'
##' @name coerce-QFeatures
##'
##' @aliases coerce,MultiAssayExperiment,QFeatures-method
##'
##' @exportMethod coerce
##'
setAs("MultiAssayExperiment", "QFeatures", function(from) {
    QFeatures(
        experiments = experiments(from),
        colData = colData(from),
        sampleMap = sampleMap(from),
        metadata = metadata(from),
        drops = from@drops,
        assayLinks = AssayLinks(names = names(from))
    )
})

##' @rdname QFeatures-class
##'
##' @exportMethod c
setMethod(
    "c", "QFeatures",
    function(x, ...) {
        ## Retrieve the assays to add
        args <- list(...)

        ## Check arguments
        if (any(sapply(args, inherits, "SummarizedExperiment")) ||
            any(sapply(args, inherits, "List")) ||
            any(sapply(args, is.list))) {
            stop(
                "Trying to combine a QFeatures object with objects that ",
                "inherit from SummarizedExperiment, List, or ",
                "list. Consider using 'addAssay()' instead."
            )
        } else if (any(sapply(args, class) == "MultiAssayExperiment")) {
            stop(
                "Trying to combine a QFeatures object with one ",
                "or more MultiAssayExperiment objects. You must ",
                "first coerce these objects to QFeatures using ",
                "'as(object, \"QFeatures\")'."
            )
        } else if (!all(sapply(args, inherits, "QFeatures"))) {
            args <- lapply(args, as, "QFeatures")
        }
        if (length(names(args))) {
            warning("Argument names are provided but will be ignored.")
        }

        ## Combine the different slots
        el <- .combineAssays(x, args)
        cd <- .combineColData(x, args)
        al <- .combineAssayLinks(x, args)

        QFeatures(
            experiments = el,
            colData = cd,
            assayLinks = al
        )
    }
)

## Internal function to combine the assays of x with the assays of each
## element in y.
## @param x A QFeatures object
## @param y A list-like object where each element is expected to be a
##     QFeatures
.combineAssays <- function(x, y) {
    Reduce(c, lapply(y, experiments), init = experiments(x))
}

## Internal function to combine the colData of x with the colData of each
## element in y.
## @param x A QFeatures object
## @param y A list-like object where each element is expected to be a
##     QFeatures
.combineColData <- function(x, y) {
    if (!length(y)) {
        return(x)
    }
    out <- colData(x)
    err <- c()
    for (i in seq_along(y)) {
        yy <- colData(y[[i]])
        cn <- .checkDataConflict(out, yy)
        if (length(cn)) {
            err <- c(err, paste0(cn, " (in argument ", i + 1, ")"))
        }
        out <- .transferData(out, yy)
    }
    if (length(err)) {
        stop(
            "Column(s) in the colData have conflicting ",
            "information when combining the QFeatures ",
            "objects. Problematic column(s): ",
            paste(err, collapse = ", ")
        )
    }
    out
}

## Internal function to combine the AssayLinks of x with the AssayLinks
## of each element in y.
## @param x A QFeatures object
## @param y A list-like object where each element is expected to be a
##     QFeatures
.combineAssayLinks <- function(x, y) {
    Reduce(c, lapply(y, attr, "assayLinks"), init = x@assayLinks)
}

##' @rdname QFeatures-class
##'
##' @param use.names A `logical(1)` indicating if the names on x
##'     should be propagated to the returned matrix or vector.
##'
##' @importFrom BiocGenerics dims
##' @exportMethod dims
setMethod(
    "dims", "QFeatures",
    function(x, use.names = TRUE) {
        vapply(experiments(x), dim, USE.NAMES = use.names, integer(2))
    }
)

##' @rdname QFeatures-class
##' @importFrom BiocGenerics nrows
##' @exportMethod nrows
setMethod(
    "nrows", "QFeatures",
    function(x, use.names = TRUE) {
        vapply(experiments(x), nrow, USE.NAMES = use.names, integer(1))
    }
)

##' @rdname QFeatures-class
##' @importFrom BiocGenerics ncols
##' @exportMethod ncols
setMethod(
    "ncols", "QFeatures",
    function(x, use.names = TRUE) {
        vapply(experiments(x), ncol, USE.NAMES = use.names, integer(1))
    }
)

##' @rdname QFeatures-class
##'
##' @param use.names A `logical(1)` indicating whether the rownames of
##'     each assay should be propagated to the corresponding `rowData`.
##'
setMethod(
    "rowData", "QFeatures",
    function(x, use.names = TRUE, ...) {
        List(lapply(experiments(x), function(xx) {
            mcols(xx, use.names = use.names, ...)
        }))
    }
)

##' @rdname QFeatures-class
##'
##' @export
setReplaceMethod(
    "rowData", c("QFeatures", "DataFrameList"),
    function(x, value) {
        i <- intersect(names(value), names(x))
        if (!length(i)) {
            warning(
                "Could not find a common assay between ",
                "'names(value)' and names(object)"
            )
            return(x)
        }
        el <- experiments(x)
        for (ii in i) {
            rowData(el[[ii]])[, colnames(value[[ii]])] <-
                value[[ii]]
        }
        BiocGenerics:::replaceSlots(x,
            ExperimentList = el,
            check = FALSE
        )
    }
)

##' @rdname QFeatures-class
##'
##' @export
setReplaceMethod(
    "rowData", c("QFeatures", "ANY"),
    function(x, value) {
        value <- endoapply(value, as, "DataFrame")
        value <- as(value, "List")
        rowData(x) <- value
        x
    }
)

##' @rdname QFeatures-class
##'
##' @export
rbindRowData <- function(object, i) {
    ## Extract the rowData and column names from the desired assay(s)
    rdlist <- rowData(object)[i]
    rdNames <- rowDataNames(object)[i]
    ## Get the common variables between the selected rowData
    commonCols <- Reduce(intersect, rdNames)
    if (!length(commonCols)) {
        warning("No common columns between rowData tables were found.")
        return(DataFrame())
    }
    ## Add assay and rowname to the rowData
    rdlist <- lapply(
        names(rdlist),
        function(x) {
            cbind(
                assay = x,
                rowname = rownames(rdlist[[x]]),
                rdlist[[x]][, commonCols]
            )
        }
    )
    ## Row bind all tables in one DataFrame
    rdlist <- do.call(rbind, rdlist)
    rownames(rdlist) <- NULL
    rdlist
}



##' @rdname QFeatures-class
##'
##' @export
selectRowData <- function(x, rowvars) {
    stopifnot(inherits(x, "QFeatures"))
    rowvars <- as.character(rowvars)
    allvars <- unique(unlist(rowDataNames(x)))
    missingvars <- setdiff(rowvars, allvars)
    if (length(missingvars)) {
        message(length(missingvars), " missing/mis-typed rowvars.")
    }
    for (i in seq_len(length(x))) {
        rd <- rowData(x[[i]])
        rowData(x[[i]]) <- rd[, colnames(rd) %in% rowvars, drop = FALSE]
    }
    x
}


##' @rdname QFeatures-class
##'
##' @importFrom Biobase fData
##'
##' @export
rowDataNames <- function(x) {
    stopifnot(inherits(x, "MultiAssayExperiment"))
    CharacterList(lapply(
        experiments(x),
        function(xx) {
            if (inherits(xx, "SummarizedExperiment")) {
                colnames(rowData(xx))
            } else if (inherits(xx, "eSet")) {
                colnames(Biobase::fData(xx))
            } else {
                NA_character_
            }
        }
    ))
}


##' @rdname QFeatures-class
##'
##' @exportMethod names<-
setReplaceMethod(
    "names", c("QFeatures", "character"),
    function(x, value) {
        key_vals <- cbind(names(x), value)
        x <- callNextMethod(x, value)
        names(x@assayLinks) <- value
        for (i in seq_len(length(x))) {
            al <- x@assayLinks[[i]]
            al@name <- unname(key_vals[key_vals[, 1] == al@name, 2])
            if (!all(is.na(al@from))) {
                al@from <- unname(key_vals[key_vals[, 1] %in% al@from, 2])
            }
            x@assayLinks[[i]] <- al
        }
        x
    }
)

##' @param y An object that inherits from `SummarizedExperiment` or a
##'     *named* list of assays. When `y` is a list, each element must
##'     inherit from a `SummarizedExperiment` and the names of the
##'     list are used as the names of the assays to add. Hence, the
##'     list names must be unique and cannot overlap with the names of
##'     the assays already present in `x`.
##'
##' @param name A `character(1)` naming the single assay. Ignored if
##'     `y` is a list of assays.
##'
##' @param assayLinks An optional [AssayLinks].
##'
##' @md
##'
##' @rdname QFeatures-class
##'
##' @export
addAssay <- function(x,
    y,
    name,
    assayLinks) {
    ## Check arguments
    stopifnot(inherits(x, "QFeatures"))
    y <- .checkAssaysToInsert(y, x, name, replace = FALSE)

    ## Check (or create) assayLinks
    if (!missing(assayLinks)) {
        if (inherits(assayLinks, "AssayLink")) {
            assayLinks <- AssayLinks(assayLinks)
        }
        if (!identical(sort(names(assayLinks)), sort(names(y)))) {
            stop("'assayLinks' must be named after the assay(s) in 'y'.")
        }
    } else {
        assayLinks <- AssayLinks(names = names(y))
    }

    ## Update the colData
    cd <- .updateColDataFromAssays(x, y)

    ## Add the assay to the ExperimentList
    ## NOTE: we replace using the `@` slot. Although not recommended,
    ## this bypasses the checks of all the elements (using
    ## `validObject`) in the ExperimentList as this is already
    ## performed when building the QFeatures object and `y` is checked
    ## at the beginning of the function. This leads to a reduction in
    ## computational time.
    el <- experiments(x)
    for (ii in names(y)) {
        el@listData[[ii]] <- y[[ii]]
    }

    ## Update the sampleMap
    smap <- MultiAssayExperiment:::.sampleMapFromData(cd, el)

    ## Update the AssayLinks
    al <- append(x@assayLinks, assayLinks)

    ## Update the QFeatures object with all the updated parts
    BiocGenerics:::replaceSlots(
        object = x,
        ExperimentList = el,
        colData = cd,
        sampleMap = smap,
        assayLinks = al,
        check = FALSE
    )
}

##' @md
##'
##' @rdname QFeatures-class
##'
##' @export
removeAssay <- function(x, i) {
    i <- .normIndex(x, i)
    x[, , !names(x) %in% i]
}

##' @md
##'
##' @rdname QFeatures-class
##'
##' @export
replaceAssay <- function(x,
    y,
    i) {
    ## Check arguments
    stopifnot(inherits(x, "QFeatures"))
    if (!missing(i)) i <- .normIndex(x, i)
    y <- .checkAssaysToInsert(y, x, i, replace = TRUE)

    ## Update the colData
    cd <- .updateColDataFromAssays(x, y)

    ## Replace the assay to the ExperimentList
    ## NOTE: we replace using the `@` slot. Although not recommended,
    ## this bypasses the checks of all the elements (using
    ## `validObject`) in the ExperimentList as this is already
    ## performed when building the QFeatures object. This leads to a
    ## reduction in computational time.
    el <- experiments(x)
    for (ii in names(y)) {
        el@listData[[ii]] <- y[[ii]]
    }

    ## Update the sampleMap
    smap <- MultiAssayExperiment:::.sampleMapFromData(cd, el)

    ## Update the AssayLinks
    al <- x@assayLinks
    allfrom <- lapply(al, function(x) x@from)
    for (ii in names(y)) {
        if (identical(
            sort(rownames(x[[ii]])),
            sort(rownames(y[[ii]]))
        ) &&
            identical(
                sort(colnames(x[[ii]])),
                sort(colnames(y[[ii]]))
            )) {
            next()
        }

        al[[ii]] <- AssayLink(ii)
        repl <- names(allfrom)[sapply(allfrom, function(x) any(x %in% ii))]
        for (jj in repl) {
            if (inherits(al[[jj]]@hits, "List")) {
                al[[jj]]@from <- al[[jj]]@from[al[[jj]]@from != ii]
                al[[jj]]@hits <- al[[jj]]@hits[names(al[[jj]]@hits) != ii]
                if (length(al[[jj]]@hits) == 1) {
                    al[[jj]]@hits <- al[[jj]]@hits[[1]]
                }
            } else {
                al[[jj]] <- AssayLink(jj)
            }
        }
    }
    if (!identical(al, x@assayLinks)) {
        warning(
            "Links between assays were lost/removed during ",
            "replacement. See '?addAssayLink' to restore them ",
            "manually. "
        )
    }

    ## Update the QFeatures object with all the updated parts
    BiocGenerics:::replaceSlots(
        object = x,
        ExperimentList = el,
        colData = cd,
        sampleMap = smap,
        assayLinks = al,
        check = FALSE
    )
}

##' @rdname QFeatures-class
##'
##' @export
setReplaceMethod(
    "[[", c("QFeatures", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., value) {
        if (length(i) != 1) {
            stop(
                "'x[[i]] <- value' does not allow multiple ",
                "replacements. Consider using 'addAssay()', ",
                "'replaceAssay()' or 'removeAssay()' instead."
            )
        }
        i <- .normIndex(x, i, allowAbsent = TRUE)
        if (!missing(j) || length(list(...))) {
            stop("invalid replacement")
        }
        if (i %in% names(x)) {
            if (is.null(value)) {
                return(removeAssay(x = x, i = i))
            } else {
                return(replaceAssay(x = x, y = value, i = i))
            }
        } else {
            return(addAssay(x = x, y = value, name = i))
        }
    }
)

## Internal function that normalize the assay indexing. In this
## context, normalization means that the returned assay index is a
## character() that complies to QFeatures assay selection.
##
## @param object A QFeatures object
##
## @param i A logical(), numeric(), factor() or character() that
##     selects an assay in object. When logical, the length of i must
##     be identical to the number of assays in object.
##
## @param allowAbsent A logical() indicating whether the i is allowed
##     to be absent from object. This argument is only applicable when
##     i is a character().
##
## @return A character() with assay names present in object, or new
##     assay names (when allowAbsent = FALSE).
.normIndex <- function(object, i, allowAbsent = FALSE) {
    if (is.logical(i) & length(i) != length(object)) {
        stop(
            "The assay index ('i') is logical but its does not ",
            "match the number of assays in the QFeatures object."
        )
    }
    if (is.factor(i)) i <- as.character(i)
    if (is.numeric(i) || is.logical(i)) {
        i <- names(object)[i]
    }
    if (!length(i)) stop("No assay selected.")
    if (any(is.na(i))) {
        stop("'i' has out of bounds entries")
    }
    if (!allowAbsent & any(mis <- !i %in% names(object))) {
        stop(
            "The following assay(s) is/are not found:",
            paste(i[mis], collapse = ",")
        )
    }
    i
}

.checkAssaysToInsert <- function(y, x, name, replace = FALSE) {
    ## Convert y to a list, if not already a list and check content
    if (!is.list(y) && !inherits(y, "List")) {
        stopifnot(is.character(name))
        y <- structure(list(y), .Names = name[1])
    } else {
        if (!missing(name)) {
            warning("'y' is provided as a list, 'name' is ignored.")
        }
        if (length(names(y)) != length(y)) {
            stop("When 'y' is a list, it must be a named List.")
        }
    }
    ## Make sure the assays comply to the requirements
    sapply(y, validObject) ## throws an error if any assay is corrupt
    if (any(duplicated(names(y)))) {
        stop("Replacement names must be unique.")
    }
    if (!replace && any(names(y) %in% names(x))) {
        stop(
            "One or more assay names are already present in 'x'. ",
            "See 'replaceAssay()' if you want to replace assays."
        )
    }
    if (replace && !all(names(y) %in% names(x))) {
        stop(
            "One or more assay names are not in 'x'. See ",
            "'addAssay()' if you want to add assays."
        )
    }
    if (!all(sapply(y, inherits, "SummarizedExperiment"))) {
        stop(
            "The replacement object(s) should inherit from ",
            "SummarizedExperiment."
        )
    }
    if (any(sapply(y, function(yy) any(duplicated(rownames(yy)))))) {
        stop("The replacement object(s) should have unique row names.")
    }
    ## Return the valid y
    y
}

## Internal function that will add rows and eventually columns in the
## colData based on a new SummarizedExperiment object
##
## @param x An instance of class [QFeatures].
## @param y A list of SummarizedExperiments containing the colData
##     that must be adapted
##
## The function returns the updated colData.
##
.updateColDataFromAssays <- function(x, y) {
    cd <- colData(x)
    ## Make sure we do not override existing colData
    err <- c()
    for (i in names(y)) {
        cn <- .checkDataConflict(cd, colData(y[[i]]))
        if (length(cn)) {
            err <- c(err, paste0(cn, " (in ", i, ")"))
        }
    }
    if (length(err) > 0) {
        stop(
            "Column(s) in the colData in y have ",
            "conflicting information with the ",
            "QFeatures colData. Problematic ",
            "column(s): ", paste(err, collapse = ", ")
        )
    }

    ## Remove lost samples (in case of replacement)
    if (any(names(y) %in% names(x))) {
        cnOld <- cnNew <- colnames(x)
        repl <- names(y)[names(y) %in% names(x)]
        for (ii in repl) {
            cnNew[[ii]] <- colnames(y[[ii]])
        }
        oldSamples <- setdiff(
            unique(unlist(cnOld)),
            unique(unlist(cnNew))
        )
        if (length(oldSamples)) {
            cd <- cd[!rownames(cd) %in% oldSamples, , drop = FALSE]
        }
    }
    ## Perform the actual colData transfer
    for (i in names(y)) {
        cd <- .transferData(cd, colData(y[[i]]))
    }
    cd
}


## Internal function that checks for data clashes between 2 tables
##
## @param x and y are data tables (DataFrame or data.frame)
##
## returns a character vector with problematic column names where a
## clash was identified. Returns an empty character vector if no problem.
.checkDataConflict <- function(x, y) {
    rn <- intersect(rownames(x), rownames(y))
    cn <- intersect(colnames(x), colnames(y))
    if (length(rn) == 0 || length(cn) == 0) {
        return(character(0))
    }
    ## We consider a problem when:
    isProbl <- sapply(cn, function(ii) {
        ## i. the overlaping colData column are different
        !identical(x[rn, ii], y[rn, ii]) &&
            ## ii. the colData x is not all missing
            !all(is.na(x[rn, ii]))
    })
    ## Return the problematic column names
    cn[isProbl]
}

## Internal function the transfers the data of y into x taking new
## rows into account
## @param x and y are data tables
##
## returns a single table where y has been transfered into x
.transferData <- function(x, y) {
    ## Add new samples names to cd and fill with NA
    newSamples <- setdiff(rownames(y), rownames(x))
    if (length(newSamples)) {
        newCd <- DataFrame(row.names = newSamples)
        newCd[, colnames(x)] <- NA
        x <- rbind(x, newCd)
    }
    ## If coldata is available, add it to cd
    if (ncol(y) != 0) {
        x[rownames(y), colnames(y)] <- y
    }
    x
}

##' @param verbose logical (default FALSE) whether to print extra messages
##'
##' @rdname QFeatures-class
##'
##' @exportMethod updateObject
setMethod(
    "updateObject", "QFeatures",
    function(object, ..., verbose = FALSE) {
        if (verbose) {
            message("updateObject(object = 'QFeatures')")
        }
        ## Update slots that are specific to QFeatures
        object@assayLinks <- updateObject(object@assayLinks, ...,
            verbose = verbose
        )
        ## Update MAE slots
        callNextMethod()
    }
)


##' @param dims `numeric()` that defines the dimensions to consider to
##'     drop empty assays. 1 for rows (i.e. assays without any
##'     features) and 2 for columns (i.e. assays without any
##'     samples). Default is `1:2`. Any value other that 1 and/or 2
##'     will trigger an error.
##'
##' @rdname QFeatures-class
##'
##' @export
dropEmptyAssays <- function(object, dims = 1:2) {
    stopifnot(inherits(object, "QFeatures"))
    if (!all(dims %in% 1:2)) {
        stop("Argument 'dims' must be in '1:2'.")
    }
    if (1 %in% dims) {
        object <- object[, , nrows(object) > 0]
    }
    if (2 %in% dims) {
        object <- object[, , ncols(object) > 0]
    }
    if (!length(object)) {
        return(QFeatures())
    }
    object
}


##' Set and Get QFeatures Type
##'
##' Developer-level functions to set and retrieve the type of a `QFeatures`
##' object. This type can help internal methods adapt their behaviour to the
##' structure of the data.
##'
##' @param object An instance of class [QFeatures].
##' @param type `character(1)` defining the type of the QFeatures.
##'     Must be one of the values returned by [validQFeaturesTypes()].
##'
##' @return
##' - `setQFeaturesType()`: returns the updated `QFeatures` object with
##'   the type stored in its metadata.
##' - `getQFeaturesType()`: returns a character string indicating the
##'   type of the `QFeatures` object. If no type is explicitly set,
##'   it is inferred from the class of the experiments.
##'   If the QFeatures contains any `SingleCellExperiment` objects,
##'   the type is set to "scp". Otherwise, it is set to "bulk".
##' - `validQFeaturesTypes()`: character vector of valid QFeatures types.
##'
##' @section Warning:
##' These functions are intended for package developers and internal use.
##' End users should typically not call them directly.
##' @details
##' These functions control an internal metadata slot (`._type`) used to
##' distinguish between different structural uses of `QFeatures` objects.
##' This slot is directly accessible with `metadata(object)[["._type"]]`.
##'
##' @note The `QFeatures` type slot was introduced because, in the
##'       context of the `scp` package, we found that `SingleCellExperiment`
##'       objects were slower than `SummarizedExperiment`
##'       objects (\href{https://github.com/UCLouvain-CBIO/scp/issues/83}{GH issue: scp#83}).
##'       As a result, we started using `SummarizedExperiment` objects
##'       within `scp`. However, to retain information about the type of
##'       data being handled, we introduced the `QFeatures` type slot.
##'       This slot is, for example, used in the `show` method of `QFeatures`.
##'
##' @rdname QFeatures-type
##' @keywords internal
##' @export
setQFeaturesType <- function(object, type) {
  stopifnot(inherits(object, "QFeatures"))
  if (!type %in% validQFeaturesTypes()) {
    stop(
      "Invalid QFeatures type. Must be one of: ",
      paste(validQFeaturesTypes(), collapse = ", ")
    )
  }
  metadata(object)[["._type"]] <- type
  object
}

##' @rdname QFeatures-type
##' @keywords internal
##' @export
getQFeaturesType <- function(object) {
  stopifnot(inherits(object, "QFeatures"))
  type <- metadata(object)[["._type"]]
  if (is.null(type)) {
    message(paste("No explicit type set for this QFeatures object,",
                  "choosing a type in fonction of experiments classes"))
    if (any(sapply(experiments(object),
                   function(x) inherits(x, "SingleCellExperiment")))) {
      type <- "scp"
    } else {
      type <- "bulk"
    }
  }
  type
}

##' @rdname QFeatures-type
##' @keywords internal
##' @export
validQFeaturesTypes <- function() {
  c("bulk", "scp")
}