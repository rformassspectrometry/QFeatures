# Reshape into a long data format

The `longForm()` method transform a
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
or SummarizedExperiment instance into a long *tidy* DataFrame that
contains the assay data, where each quantitative value is reported on a
separate line. `colData` and `rowData` varibales can also be added. This
function is an extension of the `longForm()` method in the
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html).

Note that the previous `longFormat` implementation is not defunct.

## Usage

``` r
# S4 method for class 'QFeatures'
longForm(object, colvars = NULL, rowvars = NULL, i = 1L)

# S4 method for class 'SummarizedExperiment'
longForm(object, colvars = NULL, rowvars = NULL, i = seq_along(assays(object)))
```

## Arguments

- object:

  An instance of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  or SummarizedExperiment.

- colvars:

  A [`character()`](https://rdrr.io/r/base/character.html) that selects
  column(s) in the `colData`.

- rowvars:

  A [`character()`](https://rdrr.io/r/base/character.html) with the
  names of the `rowData` variables (columns) to retain in any assay.

- i:

  When `object` is an instance of class `QFeatures`, a `numeric(1)`
  indicating what assay within each `SummarizedExperiment` object to
  return. Default is `1L`. If `object` is a `SummarizedExperiment`, a
  [`numeric()`](https://rdrr.io/r/base/numeric.html) indicating what
  assays to pull and convert. Default is to use all assays.

## Value

A `DataFrame` instance.

## Examples

``` r

data(feat2)

longForm(feat2)
#> DataFrame with 84 rows and 5 columns
#>           assay     primary     rowname     colname     value
#>     <character> <character> <character> <character> <numeric>
#> 1        assay1          S1           a          S1         1
#> 2        assay1          S1           b          S1         2
#> 3        assay1          S1           c          S1         3
#> 4        assay1          S1           d          S1         4
#> 5        assay1          S1           e          S1         5
#> ...         ...         ...         ...         ...       ...
#> 80       assay3         S12           j         S12        24
#> 81       assay3         S12           k         S12        25
#> 82       assay3         S12           l         S12        26
#> 83       assay3         S12           m         S12        27
#> 84       assay3         S12           n         S12        28

## add a colData variable and use it in longFrom
colData(feat2)$colvar <- paste0("Var", 1:12)
colData(feat2)
#> DataFrame with 12 rows and 1 column
#>          colvar
#>     <character>
#> S1         Var1
#> S2         Var2
#> S3         Var3
#> S4         Var4
#> S5         Var5
#> ...         ...
#> S8         Var8
#> S9         Var9
#> S10       Var10
#> S11       Var11
#> S12       Var12
longForm(feat2, colvars = "colvar")
#> DataFrame with 84 rows and 6 columns
#>           assay     primary     rowname     colname     value      colvar
#>     <character> <character> <character> <character> <numeric> <character>
#> 1        assay1          S1           a          S1         1        Var1
#> 2        assay1          S1           b          S1         2        Var1
#> 3        assay1          S1           c          S1         3        Var1
#> 4        assay1          S1           d          S1         4        Var1
#> 5        assay1          S1           e          S1         5        Var1
#> ...         ...         ...         ...         ...       ...         ...
#> 80       assay3         S12           j         S12        24       Var12
#> 81       assay3         S12           k         S12        25       Var12
#> 82       assay3         S12           l         S12        26       Var12
#> 83       assay3         S12           m         S12        27       Var12
#> 84       assay3         S12           n         S12        28       Var12

## use a rowData variable in longFrom
rowDataNames(feat2)
#> CharacterList of length 3
#> [["assay1"]] Prot x
#> [["assay2"]] Prot x y
#> [["assay3"]] Prot x y
longForm(feat2, rowvar = "Prot")
#> DataFrame with 84 rows and 6 columns
#>           assay     primary     rowname     colname     value        Prot
#>     <character> <character> <character> <character> <numeric> <character>
#> 1        assay1          S1           a          S1         1          Pa
#> 2        assay1          S1           b          S1         2          Pb
#> 3        assay1          S1           c          S1         3          Pc
#> 4        assay1          S1           d          S1         4          Pd
#> 5        assay1          S1           e          S1         5          Pe
#> ...         ...         ...         ...         ...       ...         ...
#> 80       assay3         S12           j         S12        24          Pj
#> 81       assay3         S12           k         S12        25          Pk
#> 82       assay3         S12           l         S12        26          Pl
#> 83       assay3         S12           m         S12        27          Pm
#> 84       assay3         S12           n         S12        28          Pn

## use both col/rowData
longForm(feat2, colvar = "colvar", rowvar = "Prot")
#> DataFrame with 84 rows and 7 columns
#>           assay     primary     rowname     colname     value      colvar
#>     <character> <character> <character> <character> <numeric> <character>
#> 1        assay1          S1           a          S1         1        Var1
#> 2        assay1          S1           b          S1         2        Var1
#> 3        assay1          S1           c          S1         3        Var1
#> 4        assay1          S1           d          S1         4        Var1
#> 5        assay1          S1           e          S1         5        Var1
#> ...         ...         ...         ...         ...       ...         ...
#> 80       assay3         S12           j         S12        24       Var12
#> 81       assay3         S12           k         S12        25       Var12
#> 82       assay3         S12           l         S12        26       Var12
#> 83       assay3         S12           m         S12        27       Var12
#> 84       assay3         S12           n         S12        28       Var12
#>            Prot
#>     <character>
#> 1            Pa
#> 2            Pb
#> 3            Pc
#> 4            Pd
#> 5            Pe
#> ...         ...
#> 80           Pj
#> 81           Pk
#> 82           Pl
#> 83           Pm
#> 84           Pn

## also works on a single SE
se <- getWithColData(feat2, 1)
#> Warning: 'experiments' dropped; see 'drops()'
longForm(se)
#> DataFrame with 40 rows and 4 columns
#>         rowname     colname     value assayName
#>     <character> <character> <numeric> <integer>
#> 1             a          S1         1         1
#> 2             b          S1         2         1
#> 3             c          S1         3         1
#> 4             d          S1         4         1
#> 5             e          S1         5         1
#> ...         ...         ...       ...       ...
#> 36            f          S4        36         1
#> 37            g          S4        37         1
#> 38            h          S4        38         1
#> 39            i          S4        39         1
#> 40            j          S4        40         1
longForm(se, colvar = "colvar")
#> DataFrame with 40 rows and 5 columns
#>         rowname     colname     value assayName      colvar
#>     <character> <character> <numeric> <integer> <character>
#> 1             a          S1         1         1        Var1
#> 2             b          S1         2         1        Var1
#> 3             c          S1         3         1        Var1
#> 4             d          S1         4         1        Var1
#> 5             e          S1         5         1        Var1
#> ...         ...         ...       ...       ...         ...
#> 36            f          S4        36         1        Var4
#> 37            g          S4        37         1        Var4
#> 38            h          S4        38         1        Var4
#> 39            i          S4        39         1        Var4
#> 40            j          S4        40         1        Var4
longForm(se, rowvar = "Prot")
#> DataFrame with 40 rows and 5 columns
#>         rowname     colname     value assayName        Prot
#>     <character> <character> <numeric> <integer> <character>
#> 1             a          S1         1         1          Pa
#> 2             b          S1         2         1          Pb
#> 3             c          S1         3         1          Pc
#> 4             d          S1         4         1          Pd
#> 5             e          S1         5         1          Pe
#> ...         ...         ...       ...       ...         ...
#> 36            f          S4        36         1          Pf
#> 37            g          S4        37         1          Pg
#> 38            h          S4        38         1          Ph
#> 39            i          S4        39         1          Pi
#> 40            j          S4        40         1          Pj
longForm(se, colvar = "colvar", rowvar = "Prot")
#> DataFrame with 40 rows and 6 columns
#>         rowname     colname     value assayName      colvar        Prot
#>     <character> <character> <numeric> <integer> <character> <character>
#> 1             a          S1         1         1        Var1          Pa
#> 2             b          S1         2         1        Var1          Pb
#> 3             c          S1         3         1        Var1          Pc
#> 4             d          S1         4         1        Var1          Pd
#> 5             e          S1         5         1        Var1          Pe
#> ...         ...         ...       ...       ...         ...         ...
#> 36            f          S4        36         1        Var4          Pf
#> 37            g          S4        37         1        Var4          Pg
#> 38            h          S4        38         1        Var4          Ph
#> 39            i          S4        39         1        Var4          Pi
#> 40            j          S4        40         1        Var4          Pj
```
