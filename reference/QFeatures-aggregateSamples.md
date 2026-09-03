# Aggregate sets' quantitative samples

This function aggregates quantitative samples, applying a summarisation
function (`fun`) to sets of samples defined by a `colData` variable.

For `QFeatures` objects, the aggregated set `colData` is added to the
global `colData` using
[`addAssay()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md).
If an aggregated group name already exists as a row name in the global
`colData`, it is reused only when shared `colData` variables do not
contain conflicting values, otherwise an error is thrown. The
aggregation of samples should take place after run joining.

## Usage

``` r
# S4 method for class 'SummarizedExperiment'
aggregateSamples(object, scol, fun = rowRobustSummary, ...)

# S4 method for class 'QFeatures'
aggregateSamples(
  object,
  i,
  scol,
  name = "newAssay",
  fun = rowRobustSummary,
  ...
)
```

## Arguments

- object:

  An instance of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  or
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

- scol:

  A [`character()`](https://rdrr.io/r/base/character.html) naming the
  `colData` variable defining how to group samples. When multiple
  variables are supplied, their values are concatenated with `_` to
  build the grouping variable.

- fun:

  A function used to aggregate the samples. The function should apply
  its summarization row-wise.

- ...:

  Additional parameters passed to `fun`.

- i:

  A `character(1)` naming the set that will be aggregated.

- name:

  A `character(1)` naming the new set.

## Value

A `QFeatures` object with an additional set or a `SummarizedExperiment`
object. The new set/`SummarizedExperiment` contains three assays the
first named `assay` containing the result of the aggregation, the second
named `aggsd` containing the standard deviation associated with each
aggregated data point and the third named `aggcounts` containing the
number of observation used for each aggregated data point.
