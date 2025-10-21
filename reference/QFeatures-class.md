# Quantitative MS QFeatures

Conceptually, a `QFeatures` object holds a set of *assays*, each
composed of a `matrix` (or `array`) containing quantitative data and row
annotations (meta-data). The number and the names of the columns
(samples) must always be the same across the assays, but the number and
the names of the rows (features) can vary. The assays are typically
defined as `SummarizedExperiment` objects. In addition, a `QFeatures`
object also uses a single `DataFrame` to annotate the samples (columns)
represented in all the matrices.

The `QFeatures` class extends the
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
and inherits all the functionality of the
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
class.

A typical use case for such `QFeatures` object is to represent
quantitative proteomics (or metabolomics) data, where different assays
represent quantitation data at the PSM (the main assay), peptide and
protein level, and where peptide values are computed from the PSM data,
and the protein-level data is calculated based on the peptide-level
values. The largest assay (the one with the highest number of features,
PSMs in the example above) is considered the main assay.

The recommended way to create `QFeatures` objects is the use the
[`readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.md)
function, that creates an instance from tabular data. The `QFeatures`
constructor can be used to create objects from their bare parts. It is
the user's responsability to make sure that these match the class
validity requirements.

## Usage

``` r
QFeatures(..., assayLinks = NULL)

# S4 method for class 'QFeatures'
show(object)

# S3 method for class 'QFeatures'
plot(x, interactive = FALSE, ...)

# S4 method for class 'QFeatures,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'QFeatures,character,ANY,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'QFeatures'
c(x, ...)

# S4 method for class 'QFeatures'
dims(x, use.names = TRUE)

# S4 method for class 'QFeatures'
nrows(x, use.names = TRUE)

# S4 method for class 'QFeatures'
ncols(x, use.names = TRUE)

# S4 method for class 'QFeatures'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'QFeatures,DataFrameList'
rowData(x) <- value

# S4 method for class 'QFeatures,ANY'
rowData(x) <- value

rbindRowData(object, i)

selectRowData(x, rowvars)

rowDataNames(x)

# S4 method for class 'QFeatures,character'
names(x) <- value

addAssay(x, y, name, assayLinks)

removeAssay(x, i)

replaceAssay(x, y, i)

# S4 method for class 'QFeatures,ANY,ANY'
x[[i, j, ...]] <- value

# S4 method for class 'QFeatures'
updateObject(object, ..., verbose = FALSE)

dropEmptyAssays(object, dims = 1:2)
```

## Arguments

- ...:

  See `MultiAssayExperiment` for details. For `plot`, further arguments
  passed to
  [`igraph::plot.igraph`](https://r.igraph.org/reference/plot.igraph.html).

- assayLinks:

  An optional
  [AssayLinks](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md).

- object:

  An instance of class QFeatures.

- x:

  An instance of class QFeatures.

- interactive:

  A `logical(1)`. If `TRUE`, an interactive graph is generated using
  `plotly`. Else, a static plot using `igraph` is generated. We
  recommend interactive exploration when the `QFeatures` object contains
  more than 50 assays.

- i:

  An indexing vector. See the corresponding section in the documentation
  for more details.

- j:

  [`character()`](https://rdrr.io/r/base/character.html),
  [`logical()`](https://rdrr.io/r/base/logical.html), or
  [`numeric()`](https://rdrr.io/r/base/numeric.html) vector for
  subsetting by `colData` rows.

- drop:

  logical (default `TRUE`) whether to drop empty assay elements in the
  `ExperimentList`.

- k:

  [`character()`](https://rdrr.io/r/base/character.html),
  [`logical()`](https://rdrr.io/r/base/logical.html), or
  [`numeric()`](https://rdrr.io/r/base/numeric.html) vector for
  subsetting by assays

- use.names:

  A `logical(1)` indicating whether the rownames of each assay should be
  propagated to the corresponding `rowData`.

- value:

  The values to use as a replacement. See the corresponding section in
  the documentation for more details.

- rowvars:

  A [`character()`](https://rdrr.io/r/base/character.html) with the
  names of the `rowData` variables (columns) to retain in any assay.

- y:

  An object that inherits from `SummarizedExperiment` or a *named* list
  of assays. When `y` is a list, each element must inherit from a
  `SummarizedExperiment` and the names of the list are used as the names
  of the assays to add. Hence, the list names must be unique and cannot
  overlap with the names of the assays already present in `x`.

- name:

  A `character(1)` naming the single assay. Ignored if `y` is a list of
  assays.

- verbose:

  logical (default FALSE) whether to print extra messages

- dims:

  [`numeric()`](https://rdrr.io/r/base/numeric.html) that defines the
  dimensions to consider to drop empty assays. 1 for rows (i.e. assays
  without any features) and 2 for columns (i.e. assays without any
  samples). Default is `1:2`. Any value other that 1 and/or 2 will
  trigger an error.

## Constructors

- `QFeatures(..., assayLinks)` allows the manual construction of
  objects. It is the user's responsability to make sure these comply.
  The arguments in `...` are those documented in
  [`MultiAssayExperiment::MultiAssayExperiment()`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html).
  For details about `assayLinks`, see
  [AssayLinks](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md).
  An example is shown below.

- The
  [`readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.md)
  function constructs a `QFeatures` object from text-based spreadsheet
  or a `data.frame` used to generate an assay. See the function manual
  page for details and an example.

## Accessors

- The `QFeatures` class extends the
  [MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
  class and inherits all its accessors and replacement methods.

- The `rowData` method returns a `DataFrameList` containing the
  `rowData` for each assay of the `QFeatures` object. On the other hand,
  `rowData` can be modified using `rowData(x) <- value`, where `value`
  is a list of tables that can be coerced to `DFrame` tables. The names
  of `value` point to the assays for which the `rowData` must be
  replaced. The column names of each table are used to replace the data
  in the existing `rowData`. If the column name does not exist, a new
  column is added to the `rowData`.

- The `rbindRowData` functions returns a `DFrame` table that contains
  the row binded `rowData` tables from the selected assays. In this
  context, `i` is a
  [`character()`](https://rdrr.io/r/base/character.html),
  [`integer()`](https://rdrr.io/r/base/integer.html) or
  [`logical()`](https://rdrr.io/r/base/logical.html) object for
  subsetting assays. Only rowData variables that are common to all
  assays are kept.

- The `rowDataNames` accessor returns a list with the `rowData` variable
  names.

- The
  [`longForm()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-longForm.md)
  accessor takes a `QFeatures` instance and returns it in a long *tidy*
  `DataFrame`, where each quantitative value is reported on a separate
  line.

## Adding, removing and replacing assays

- The
  [`aggregateFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  function creates a new assay by aggregating features of an existing
  assay.

- `addAssay(x, y, name, assayLinks)`: Adds one or more new assay(s) `y`
  to the `QFeatures` instance `x`. `name` is a `character(1)` naming the
  assay if only one assay is provided, and is ignored if `y` is a list
  of assays. `assayLinks` is an optional
  [AssayLinks](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md).
  The `colData(y)` is automatically added to `colData(x)` by matching
  sample names, that is `colnames(y)`. If the samples are not present in
  `x`, the rows of `colData(x)` are extended to account for the new
  samples. Be aware that conflicting information between the
  `colData(y)` and the `colData(x)` will result in an error.

- `removeAssay(x, i)`: Removes one or more assay(s) from the `QFeatures`
  instance `x`. In this context, `i` is a
  [`character()`](https://rdrr.io/r/base/character.html),
  [`integer()`](https://rdrr.io/r/base/integer.html) or
  [`logical()`](https://rdrr.io/r/base/logical.html) that indicates
  which assay(s) to remove.

- `replaceAssay(x, y, i)`: Replaces one or more assay(s) from the
  `QFeatures` instance `x`. In this context, `i` is a
  [`character()`](https://rdrr.io/r/base/character.html),
  [`integer()`](https://rdrr.io/r/base/integer.html) or
  [`logical()`](https://rdrr.io/r/base/logical.html) that indicates
  which assay(s) to replace. The `AssayLinks` from or to any replaced
  assays are automatically removed, unless the replacement has the same
  dimension names (columns and row, order agnostic). Be aware that
  conflicting information between `colData(y)` and `colData(x)` will
  result in an error.

- `x[[i]] <- value`: a generic method for adding (when `i` is not in
  `names(x)`), removing (when `value` is null) or replacing (when `i` is
  in `names(x)`). Note that the arguments `j` and `...` from the S4
  replacement method signature are not allowed.

## Subsetting

- QFeatures object can be subset using the `x[i, j, k, drop = TRUE]`
  paradigm. In this context, `i` is a
  [`character()`](https://rdrr.io/r/base/character.html),
  [`integer()`](https://rdrr.io/r/base/integer.html),
  [`logical()`](https://rdrr.io/r/base/logical.html) or `GRanges()`
  object for subsetting by rows. See the argument descriptions for
  details on the remaining arguments.

- The
  [`subsetByFeature()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-subsetBy.md)
  function can be used to subset a `QFeatures` object using one or
  multiple feature names that will be matched across different assays,
  taking the aggregation relation between assays.

- The `selectRowData(x, rowvars)` function can be used to select a
  limited number of `rowData` columns of interest named in `rowvars` in
  the `x` instance of class `QFeatures`. All other variables than
  `rowvars` will be dropped. In case an element in `rowvars` isn't found
  in any `rowData` variable, a message is printed.

- The `dropEmptyAssays(object, dims)` function removes empty assays from
  a `QFeatures`. Empty assays are defined as having 0 rows and/or 0
  columns, as defined by the `dims` argument.

## See also

- The
  [`readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.md)
  constructor and the
  [`aggregateFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  function. The *QFeatures* vignette provides an extended example.

- The
  [QFeatures-filtering](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-filtering.md)
  manual page demonstrates how to filter features based on their
  rowData.

- The
  [missing-data](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  manual page to manage missing values in `QFeatures` objects.

- The
  [QFeatures-processing](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  and
  [`aggregateFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  manual pages and *Processing* vignette describe common quantitative
  data processing methods using in quantitative proteomics.

## Author

Laurent Gatto

## Examples

``` r
## ------------------------
## An empty QFeatures object
## ------------------------

QFeatures()
#> An empty instance of class QFeatures (type: bulk)

## -----------------------------------
## Creating a QFeatures object manually
## -----------------------------------

## two assays (matrices) with matching column names
m1 <- matrix(1:40, ncol = 4)
m2 <- matrix(1:16, ncol = 4)
sample_names <- paste0("S", 1:4)
colnames(m1) <- colnames(m2) <- sample_names
rownames(m1) <- letters[1:10]
rownames(m2) <- letters[1:4]

## two corresponding feature metadata with appropriate row names
df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
                 row.names = rownames(m1))
df2 <- DataFrame(row.names = rownames(m2))

(se1 <- SummarizedExperiment(m1, df1))
#> class: SummarizedExperiment 
#> dim: 10 4 
#> metadata(0):
#> assays(1): ''
#> rownames(10): a b ... i j
#> rowData names(2): Fa Fb
#> colnames(4): S1 S2 S3 S4
#> colData names(0):
(se2 <- SummarizedExperiment(m2, df2))
#> class: SummarizedExperiment 
#> dim: 4 4 
#> metadata(0):
#> assays(1): ''
#> rownames(4): a b c d
#> rowData names(0):
#> colnames(4): S1 S2 S3 S4
#> colData names(0):

## Sample annotation (colData)
cd <- DataFrame(Var1 = rnorm(4),
                Var2 = LETTERS[1:4],
                row.names = sample_names)

el <- list(assay1 = se1, assay2 = se2)
fts1 <- QFeatures(el, colData = cd)
fts1
#> An instance of class QFeatures (type: bulk) with 2 sets:
#> 
#>  [1] assay1: SummarizedExperiment with 10 rows and 4 columns 
#>  [2] assay2: SummarizedExperiment with 4 rows and 4 columns 
fts1[[1]]
#> class: SummarizedExperiment 
#> dim: 10 4 
#> metadata(0):
#> assays(1): ''
#> rownames(10): a b ... i j
#> rowData names(2): Fa Fb
#> colnames(4): S1 S2 S3 S4
#> colData names(0):
fts1[["assay1"]]
#> class: SummarizedExperiment 
#> dim: 10 4 
#> metadata(0):
#> assays(1): ''
#> rownames(10): a b ... i j
#> rowData names(2): Fa Fb
#> colnames(4): S1 S2 S3 S4
#> colData names(0):

## Rename assay
names(fts1) <- c("se1", "se2")

## Add an assay
fts1 <- addAssay(fts1, se1[1:2, ], name = "se3")

## Get the assays feature metadata
rowData(fts1)
#> DataFrameList of length 3
#> names(3): se1 se2 se3

## Keep only the Fa variable
selectRowData(fts1, rowvars = "Fa")
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] se1: SummarizedExperiment with 10 rows and 4 columns 
#>  [2] se2: SummarizedExperiment with 4 rows and 4 columns 
#>  [3] se3: SummarizedExperiment with 2 rows and 4 columns 

## -----------------------------------
## See ?readQFeatures to create a
## QFeatures object from a data.frame
## or spreadsheet.
## -----------------------------------
```
