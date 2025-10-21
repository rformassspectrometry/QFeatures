# QFeatures processing

This manual page describes common quantitative proteomics data
processing methods using
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
objects. In the following functions, if `object` is of class
`QFeatures`, and optional assay index or name `i` can be specified to
define the assay (by name of index) on which to operate.

The following functions are currently available:

- `logTransform(object, base = 2, i, pc = 0)` log-transforms (with an
  optional pseudocount offset) the assay(s).

- `normalize(object, method, i)` normalises the assay(s) according to
  `method` (see Details).

- `scaleTransform(object, center = TRUE, scale = TRUE, i)` applies
  [`base::scale()`](https://rdrr.io/r/base/scale.html) to
  `SummarizedExperiment` and `QFeatures` objects.

- `sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)` sweeps
  out array summaries from `SummarizedExperiment` and `QFeatures`
  objects. See [`base::sweep()`](https://rdrr.io/r/base/sweep.html) for
  details.

See the *Processing* vignette for examples.

## Usage

``` r
# S4 method for class 'SummarizedExperiment'
logTransform(object, base = 2, pc = 0)

# S4 method for class 'QFeatures'
logTransform(object, i, name = "logAssay", base = 2, pc = 0)

# S4 method for class 'SummarizedExperiment'
scaleTransform(object, center = TRUE, scale = TRUE)

# S4 method for class 'QFeatures'
scaleTransform(object, i, name = "scaledAssay", center = TRUE, scale = TRUE)

# S4 method for class 'SummarizedExperiment'
normalize(object, method, ...)

# S4 method for class 'QFeatures'
normalize(object, i, name = "normAssay", method, ...)

# S4 method for class 'SummarizedExperiment'
sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)

# S4 method for class 'QFeatures'
sweep(
  x,
  MARGIN,
  STATS,
  FUN = "-",
  check.margin = TRUE,
  ...,
  i,
  name = "sweptAssay"
)
```

## Arguments

- object:

  An object of class `QFeatures` or `SummarizedExperiment`.

- base:

  `numeric(1)` providing the base with respect to which logarithms are
  computed. Defaults is 2.

- pc:

  `numeric(1)` with a pseudocount to add to the quantitative data.
  Useful when (true) 0 are present in the data. Default is 0 (no
  effect).

- i:

  A numeric vector or a character vector giving the index or the name,
  respectively, of the assay(s) to be processed.

- name:

  A `character(1)` naming the new assay name. Defaults are `logAssay`
  for `logTransform`, `scaledAssay` for `scaleTranform` and `normAssay`
  for `normalize`.

- center:

  `logical(1)` (default is `TRUE`) value or numeric-alike vector of
  length equal to the number of columns of `object`. See
  [`base::scale()`](https://rdrr.io/r/base/scale.html) for details.

- scale:

  `logical(1)` (default is `TRUE`) or a numeric-alike vector of length
  equal to the number of columns of `object`. See
  [`base::scale()`](https://rdrr.io/r/base/scale.html) for details.

- method:

  `character(1)` defining the normalisation method to apply. See
  Details.

- ...:

  Additional parameters passed to inner functions.

- x:

  An object of class `QFeatures` or `SummarizedExperiment` in `sweep`.

- MARGIN:

  As in [`base::sweep()`](https://rdrr.io/r/base/sweep.html), a vector
  of indices giving the extent(s) of `x` which correspond to `STATS`.

- STATS:

  As in [`base::sweep()`](https://rdrr.io/r/base/sweep.html), the
  summary statistic which is to be swept out.

- FUN:

  As in [`base::sweep()`](https://rdrr.io/r/base/sweep.html), the
  function to be used to carry out the sweep.

- check.margin:

  As in [`base::sweep()`](https://rdrr.io/r/base/sweep.html), a
  `logical`. If `TRUE` (the default), warn if the length or dimensions
  of `STATS` do not match the specified dimensions of `x`. Set to
  `FALSE` for a small speed gain when you know that dimensions match.

## Value

An processed object of the same class as `x` or `object`.

## Details

The `method` parameter in `normalize` can be one of `"sum"`, `"max"`,
`"center.mean"`, `"center.median"`, `"div.mean"`, `"div.median"`,
`"diff.median"`, `"quantiles`", `"quantiles.robust`" or `"vsn"`. The
[`MsCoreUtils::normalizeMethods()`](https://rdrr.io/pkg/MsCoreUtils/man/normalize.html)
function returns a vector of available normalisation methods.

- For `"sum"` and `"max"`, each feature's intensity is divided by the
  maximum or the sum of the feature respectively. These two methods are
  applied along the features (rows).

- `"center.mean"` and `"center.median"` center the respective sample
  (column) intensities by subtracting the respective column means or
  medians. `"div.mean"` and `"div.median"` divide by the column means or
  medians. These are equivalent to `sweep`ing the column means (medians)
  along `MARGIN = 2` with `FUN = "-"` (for `"center.*"`) or `FUN = "/"`
  (for `"div.*"`).

- `"diff.median"` centers all samples (columns) so that they all match
  the grand median by subtracting the respective columns medians
  differences to the grand median.

- Using `"quantiles"` or `"quantiles.robust"` applies (robust) quantile
  normalisation, as implemented in
  [`preprocessCore::normalize.quantiles()`](https://rdrr.io/pkg/preprocessCore/man/normalize.quantiles.html)
  and
  [`preprocessCore::normalize.quantiles.robust()`](https://rdrr.io/pkg/preprocessCore/man/normalize.quantiles.robust.html).
  `"vsn"` uses the
  [`vsn::vsn2()`](https://rdrr.io/pkg/vsn/man/vsn2.html) function. Note
  that the latter also glog-transforms the intensities. See respective
  manuals for more details and function arguments.

For further details and examples about normalisation, see
[`MsCoreUtils::normalize_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/normalize.html).

## Examples

``` r

MsCoreUtils::normalizeMethods()
#>  [1] "sum"              "max"              "center.mean"      "center.median"   
#>  [5] "div.mean"         "div.median"       "diff.median"      "quantiles"       
#>  [9] "quantiles.robust" "vsn"             
```
