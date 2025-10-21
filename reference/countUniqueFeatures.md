# Count Unique Features

This function counts the number of unique features per sample. A
grouping structure can be provided to count higher level features from
assays, for example counting the number of unique proteins from PSM
data.

## Usage

``` r
countUniqueFeatures(object, i, groupBy = NULL, colDataName = "count")
```

## Arguments

- object:

  An object of class `QFeatures`.

- i:

  A [`numeric()`](https://rdrr.io/r/base/numeric.html) or
  [`character()`](https://rdrr.io/r/base/character.html) vector
  indicating from which assays the `rowData` should be taken.

- groupBy:

  A `character(1)` indicating the variable name in the `rowData` that
  contains the grouping variable, for instance to count the unique
  number of peptides or proteins expressed in each samples (column). If
  `groupBy` is missing, the number of non zero elements per sample will
  be stored.

- colDataName:

  A `character(1)` giving the name of the new variable in the `colData`
  where the number of unique features will be stored. The name cannot
  already exist in the `colData`.

## Value

An object of class `QFeatures`.

## Examples

``` r
data("ft_na")
## Count number of (non-missing) PSMs
ft_na <- countUniqueFeatures(ft_na, 
                             i = "na", 
                             colDataName = "counts")
ft_na$counts
#> [1] 2 3 4
## Count number of unique rowData feature
ft_na <- countUniqueFeatures(ft_na, 
                             i = "na", 
                             groupBy = "Y",
                             colDataName = "Y_counts")
ft_na$Y_counts
#> [1] 2 2 2
```
