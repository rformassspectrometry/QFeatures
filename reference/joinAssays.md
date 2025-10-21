# Join assays in a QFeatures object

This function applies a full-join type of operation on 2 or more assays
in a `QFeatures` instance.

## Usage

``` r
joinAssays(x, i, name = "joinedAssay")
```

## Arguments

- x:

  An instance of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md).

- i:

  The indices or names of al least two assays to be joined.

- name:

  A `character(1)` naming the new assay. Default is `joinedAssay`. Note
  that the function will fail if there's already an assay with `name`.

## Value

A `QFeatures` object with an additional assay.

## Details

The rows to be joined are chosen based on the rownames of the respective
assays. It is the user's responsability to make sure these are
meaningful, such as for example refering to unique peptide sequences or
proteins.

The join operation acts along the rows and expects the samples (columns)
of the assays to be disjoint, i.e. the assays mustn't share any samples.
Rows that aren't present in an assay are set to `NA` when merged.

The `rowData` slots are also joined. However, only columns that are
shared and that have the same values for matching columns/rows are
retained. For example of a feature variable `A` in sample `S1` contains
value `a1` and variable `A` in sample `S2` in a different assay contains
`a2`, then the feature variable `A` is dropped in the merged assay.

The joined assay is linked to its parent assays through an `AssayLink`
object. The link between the child assay and the parent assays is based
on the assay row names, just like the procedure for joining the parent
assays.

## Author

Laurent Gatto

## Examples

``` r

## -----------------------------------------------
## An example QFeatures with 3 assays to be joined
## -----------------------------------------------
data(feat2)
feat2
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] assay1: SummarizedExperiment with 10 rows and 4 columns 
#>  [2] assay2: SummarizedExperiment with 4 rows and 4 columns 
#>  [3] assay3: SummarizedExperiment with 7 rows and 4 columns 

feat2 <- joinAssays(feat2, 1:3)

## Individual assays to be joined, each with 4 samples and a
## variable number of rows.
assay(feat2[[1]])
#>   S1 S2 S3 S4
#> a  1 11 21 31
#> b  2 12 22 32
#> c  3 13 23 33
#> d  4 14 24 34
#> e  5 15 25 35
#> f  6 16 26 36
#> g  7 17 27 37
#> h  8 18 28 38
#> i  9 19 29 39
#> j 10 20 30 40
assay(feat2[[2]])
#>   S5 S6 S7 S8
#> h  1  5  9 13
#> i  2  6 10 14
#> j  3  7 11 15
#> k  4  8 12 16
assay(feat2[[3]])
#>   S9 S10 S11 S12
#> a  1   8  15  22
#> b  2   9  16  23
#> j  3  10  17  24
#> k  4  11  18  25
#> l  5  12  19  26
#> m  6  13  20  27
#> n  7  14  21  28

## The joined assay contains 14 rows (corresponding to the union
## of those in the initial assays) and 12 samples
assay(feat2[["joinedAssay"]])
#>   S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12
#> j 10 20 30 40  3  7 11 15  3  10  17  24
#> a  1 11 21 31 NA NA NA NA  1   8  15  22
#> b  2 12 22 32 NA NA NA NA  2   9  16  23
#> k NA NA NA NA  4  8 12 16  4  11  18  25
#> h  8 18 28 38  1  5  9 13 NA  NA  NA  NA
#> i  9 19 29 39  2  6 10 14 NA  NA  NA  NA
#> d  4 14 24 34 NA NA NA NA NA  NA  NA  NA
#> e  5 15 25 35 NA NA NA NA NA  NA  NA  NA
#> f  6 16 26 36 NA NA NA NA NA  NA  NA  NA
#> c  3 13 23 33 NA NA NA NA NA  NA  NA  NA
#> g  7 17 27 37 NA NA NA NA NA  NA  NA  NA
#> l NA NA NA NA NA NA NA NA  5  12  19  26
#> m NA NA NA NA NA NA NA NA  6  13  20  27
#> n NA NA NA NA NA NA NA NA  7  14  21  28

## The individual rowData to be joined.
rowData(feat2[[1]])
#> DataFrame with 10 rows and 2 columns
#>          Prot         x
#>   <character> <numeric>
#> a          Pa  2.067646
#> b          Pb -0.188981
#> c          Pc  0.266870
#> d          Pd  1.671331
#> e          Pe -1.857170
#> f          Pf  1.166811
#> g          Pg  0.316521
#> h          Ph  0.976154
#> i          Pi  0.117673
#> j          Pj  2.420734
rowData(feat2[[2]])
#> DataFrame with 4 rows and 3 columns
#>          Prot          x         y
#>   <character>  <numeric> <numeric>
#> h          Ph -0.0463617  1.547125
#> i          Pi  0.6833728 -1.513248
#> j          Pj  0.7917619  0.574492
#> k          Pk -0.5226943  1.404666
rowData(feat2[[3]])
#> DataFrame with 7 rows and 3 columns
#>          Prot          x         y
#>   <character>  <numeric> <numeric>
#> a          Pa -0.0762406 -0.420478
#> b          Pb  0.5580027 -0.389534
#> j          Pj  2.5905327 -0.377114
#> k          Pk  0.1605870 -0.418806
#> l          Pl  0.7486858  1.967733
#> m          Pm  1.2365642 -0.516798
#> n          Pn -0.4848651  0.263387

## Only the 'Prot' variable is retained because it is shared among
## all assays and the values and coherent across samples (the
## value of 'Prot' for row 'j' is always 'Pj'). The variable 'y' is
## missing in 'assay1' and while variable 'x' is present is all
## assays, the values for the shared rows are different.
rowData(feat2[["joinedAssay"]])
#> DataFrame with 14 rows and 1 column
#>            Prot
#>     <character>
#> j            Pj
#> a            Pa
#> b            Pb
#> k            Pk
#> h            Ph
#> ...         ...
#> c            Pc
#> g            Pg
#> l            Pl
#> m            Pm
#> n            Pn
```
