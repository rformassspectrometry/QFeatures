# Example `QFeatures`

`feat3` is a small `QFeatures` object that contains 2 PSM-level assays
used to illustrate to creation of unique precursor identifiers and
merging, as shown in
[`createPrecursorId()`](https://rformassspectrometry.github.io/QFeatures/reference/createPrecursorId.md).

## Usage

``` r
feat4
```

## Format

An object of class `QFeatures` of length 2.

## Source

`feat4` was built from `feat3`. The source code is available in
[`inst/scripts/make-feat4.R`](https://github.com/rformassspectrometry/QFeatures/blob/master/inst/scripts/make-feat4.R)

## See also

See
[`?feat1`](https://rformassspectrometry.github.io/QFeatures/reference/feat1.md)
for other example/test data sets.

## Examples

``` r

data("feat4")
feat4
#> An instance of class QFeatures (type: bulk) with 2 sets:
#> 
#>  [1] PSM1: SummarizedExperiment with 7 rows and 2 columns 
#>  [2] PSM2: SummarizedExperiment with 8 rows and 2 columns 
```
