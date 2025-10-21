# Managing missing data

This manual page describes the handling of missing values in
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
objects. In the following functions, if `object` is of class
`QFeatures`, an optional assay index or name `i` can be specified to
define the assay (by name of index) on which to operate.

The following functions are currently available:

- `zeroIsNA(object, i)` replaces all 0 in `object` by `NA`. This is
  often necessary when third-party software assume that features that
  weren't quantified should be assigned an intensity of 0.

- `infIsNA(object, i)` replaces all infinite values in `object` by `NA`.
  This is necessary when third-party software divide expression data by
  zero values, for instance during custom normalization.

- `nNA(object, i)` returns a list of missing value summaries. The first
  element `nNA` gives a `DataFrame` with the number and the proportion
  of missing values for the whole assay; the second element `nNArows`
  provides a `DataFrame` with the number and the proportion of missing
  values for the features (rows) of the assay(s); the third element
  `nNAcols` provides the number and the proportions of missing values in
  each sample of the assay(s). When `object` has class `QFeatures` and
  additional column with the assays is provided in each element's
  `DataFrame`.

- `filterNA(object, pNA, i)` removes features (rows) that contain a
  proportion of more missing values of `pNA` or higher.

See the *Processing* vignette for examples.

## Usage

``` r
# S4 method for class 'SummarizedExperiment,missing'
zeroIsNA(object, i)

# S4 method for class 'QFeatures,integer'
zeroIsNA(object, i)

# S4 method for class 'QFeatures,numeric'
zeroIsNA(object, i)

# S4 method for class 'QFeatures,character'
zeroIsNA(object, i)

# S4 method for class 'SummarizedExperiment,missing'
infIsNA(object, i)

# S4 method for class 'QFeatures,integer'
infIsNA(object, i)

# S4 method for class 'QFeatures,numeric'
infIsNA(object, i)

# S4 method for class 'QFeatures,character'
infIsNA(object, i)

# S4 method for class 'SummarizedExperiment,missing'
nNA(object, i)

# S4 method for class 'QFeatures,integer'
nNA(object, i)

# S4 method for class 'QFeatures,numeric'
nNA(object, i)

# S4 method for class 'QFeatures,character'
nNA(object, i)

# S4 method for class 'SummarizedExperiment'
filterNA(object, pNA = 0)

# S4 method for class 'QFeatures'
filterNA(object, pNA = 0, i)
```

## Arguments

- object:

  An object of class `QFeatures` or `SummarizedExperiment`.

- i:

  One or more indices or names of the assay(s) to be processed.

- pNA:

  `numeric(1)` providing the maximum proportion of missing values per
  feature (row) that is acceptable. Feature with higher proportions are
  removed. If 0 (default), features that contain any number of `NA`
  values are dropped.

## Value

An instance of the same class as `object`.

## See also

The
[`impute()`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md)
for `QFeautres` instances.

## Examples

``` r
data(ft_na)

## Summary if missing values
nNA(ft_na, 1)
#> $nNA
#> DataFrame with 1 row and 3 columns
#>         assay       nNA       pNA
#>   <character> <integer> <numeric>
#> 1          na         3      0.25
#> 
#> $nNArows
#> DataFrame with 4 rows and 4 columns
#>         assay        name       nNA       pNA
#>   <character> <character> <integer> <numeric>
#> 1          na           a         1  0.333333
#> 2          na           b         0  0.000000
#> 3          na           c         1  0.333333
#> 4          na           d         1  0.333333
#> 
#> $nNAcols
#> DataFrame with 3 rows and 4 columns
#>         assay        name       nNA       pNA
#>   <character> <character> <integer> <numeric>
#> 1          na           A         2      0.50
#> 2          na           B         1      0.25
#> 3          na           C         0      0.00
#> 

## Remove rows with missing values
assay(filterNA(ft_na, i = 1))
#>   A B  C
#> b 2 6 10

## Replace NAs by zero and back
ft_na <- impute(ft_na, i = 1, method = "zero")
assay(ft_na)
#>    A  B  C
#> a NA  5  9
#> b  2  6 10
#> c  3 NA 11
#> d NA  8 12
ft_na <- zeroIsNA(ft_na, 1)
assay(ft_na)
#>    A  B  C
#> a NA  5  9
#> b  2  6 10
#> c  3 NA 11
#> d NA  8 12
```
