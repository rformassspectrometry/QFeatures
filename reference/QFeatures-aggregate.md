# Aggregate assays' quantitative features

This function aggregates the quantitative features of one or multiple
assays, applying a summarisation function (`fun`) to sets of features.
The `fcol` variable name points to a rowData column that defines how to
group the features during aggregate. This variable can eigher be a
vector (we then refer to an *aggregation by vector*) or an adjacency
matrix (*aggregation by matrix*).

The rowData of the aggregated `SummarizedExperiment` assays contains a
`.n` variable that provides the number of parent features that were
aggregated.

When aggregating with a vector, the newly aggregated
`SummarizedExperiment` assays also contains a new `aggcounts` assay
containing the aggregation counts matrix, i.e. the number of features
that were aggregated for each sample, which can be accessed with the
`aggcounts()` accessor.

Only the rowData columns that are invariant within a group across all
assays will be retained in the new assays' rowData.

## Usage

``` r
# S4 method for class 'QFeatures'
aggregateFeatures(
  object,
  i,
  fcol,
  name = "newAssay",
  fun = MsCoreUtils::robustSummary,
  ...
)

# S4 method for class 'SummarizedExperiment'
aggregateFeatures(object, fcol, fun = MsCoreUtils::robustSummary, ...)

# S4 method for class 'QFeatures'
adjacencyMatrix(object, i, adjName = "adjacencyMatrix")

adjacencyMatrix(object, i, adjName = "adjacencyMatrix") <- value

# S4 method for class 'SummarizedExperiment'
aggcounts(object, ...)
```

## Arguments

- object:

  An instance of class `SummarizedExperiment` or `QFeatures`.

- i:

  When adding an adjacency matrix to an assay of a `QFeatures` object,
  the index or name of the assay the adjacency matrix will be added to.
  Ignored when `x` is an `SummarizedExperiment`.

- fcol:

  A `character(1)` naming a rowdata variable (of assay `i` in case of a
  `QFeatures`) defining how to aggregate the features of the assays.
  This variable is either a `character` or a (possibly sparse) matrix.
  See below for details.

- name:

  A [`character()`](https://rdrr.io/r/base/character.html) naming the
  new assays. `name` must have the same length as i. Default is
  `newAssay`. Note that the function will fail if there's already an
  assay with `name`.

- fun:

  A function used for quantitative feature aggregation. See Details for
  examples.

- ...:

  Additional parameters passed the `fun`.

- adjName:

  `character(1)` with the variable name containing the adjacency matrix.
  Default is `"adjacencyMatrix"`.

- value:

  An adjacency matrix with row and column names. The matrix will be
  coerced to compressed, column-oriented sparse matrix (class
  `dgCMatrix`) as defined in the `Matrix` package, as generaled by the
  `sparseMatrix()` constructor.

## Value

A `QFeatures` object with an additional assay or a
`SummarizedExperiment` object (or subclass thereof).

## Details

Aggregation is performed by a function that takes a matrix as input and
returns a vector of length equal to `ncol(x)`. Examples thereof are

- [`MsCoreUtils::medianPolish()`](https://rdrr.io/pkg/MsCoreUtils/man/medianPolish.html)
  to fits an additive model (two way decomposition) using Tukey's median
  polish\_ procedure using
  [`stats::medpolish()`](https://rdrr.io/r/stats/medpolish.html);

- [`MsCoreUtils::robustSummary()`](https://rdrr.io/pkg/MsCoreUtils/man/robustSummary.html)
  to calculate a robust aggregation using
  [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) (default);

- [`base::colMeans()`](https://rdrr.io/r/base/colSums.html) to use the
  mean of each column;

- `colMeansMat(x, MAT)` to aggregate feature by the calculating the mean
  of peptide intensities via an adjacency matrix. Shared peptides are
  re-used multiple times.

- [`matrixStats::colMedians()`](https://rdrr.io/pkg/matrixStats/man/rowMedians.html)
  to use the median of each column.

- [`base::colSums()`](https://rdrr.io/r/base/colSums.html) to use the
  sum of each column;

- `colSumsMat(x, MAT)` to aggregate feature by the summing the peptide
  intensities for each protein via an adjacency matrix. Shared peptides
  are re-used multiple times.

See
[`MsCoreUtils::aggregate_by_vector()`](https://rdrr.io/pkg/MsCoreUtils/man/aggregate.html)
for more aggregation functions.

## Missing quantitative values

Missing quantitative values have different effects based on the
aggregation method employed:

- The aggregation functions should be able to deal with missing values
  by either ignoring or propagating them. This is often done with an
  `na.rm` argument, that can be passed with `...`. For example,
  `rowSums`, `rowMeans`, `rowMedians`, ... will ignore `NA` values with
  `na.rm = TRUE`, as illustrated below.

- Missing values will result in an error when using `medpolish`, unless
  `na.rm = TRUE` is used. Note that this option relies on implicit
  assumptions and/or performes an implicit imputation: when summing, the
  values are implicitly imputed by 0, assuming that the `NA` represent a
  trully absent features; when averaging, the assumption is that the
  `NA` represented a genuinely missing value.

- When using robust summarisation, individual missing values are
  excluded prior to fitting the linear model by robust regression. To
  remove all values in the feature containing the missing values, use
  [`filterNA()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md).

More generally, missing values often need dedicated handling such as
filtering (see
[`filterNA()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md))
or imputation (see
[`impute()`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md)).

## Missing values in the row data

Missing values in the row data of an assay will also impact the
resulting (aggregated) assay row data, as illustrated in the example
below. Any feature variables (a column in the row data) containing `NA`
values will be dropped from the aggregated row data. The reasons
underlying this drop are detailed in the
[`reduceDataFrame()`](https://rformassspectrometry.github.io/QFeatures/reference/reduceDataFrame.md)
manual page: only invariant aggregated rows, i.e. rows resulting from
the aggregation from identical variables, are preserved during
aggregations.

The situation illustrated below should however only happen in rare cases
and should often be imputable using the value of the other aggregation
rows before aggregation to preserve the invariant nature of that column.
In cases where an `NA` is present in an otherwise variant column, the
column would be dropped anyway.

## Using an adjacency matrix

When considering non-unique peptides explicitly, i.e. peptides that map
to multiple proteins rather than as a protein group, it is convenient to
encode this ambiguity explicitly using a peptide-by-proteins (sparse)
adjacency matrix. This matrix is typically stored in the rowdata and
set/retrieved with the `adjacencyMatrix()` function. It can be created
manually (as illustrated below) or using
`PSMatch::makeAdjacencyMatrix()`.

## See also

The *QFeatures* vignette provides an extended example and the
*Processing* vignette, for a complete quantitative proteomics data
processing pipeline. The
[`MsCoreUtils::aggregate_by_vector()`](https://rdrr.io/pkg/MsCoreUtils/man/aggregate.html)
manual page provides further details.

## Examples

``` r

## ---------------------------------------
## An example QFeatures with PSM-level data
## ---------------------------------------
data(feat1)
feat1
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 

## Aggregate PSMs into peptides
feat1 <- aggregateFeatures(feat1, "psms", "Sequence", name = "peptides")
#> Aggregated: 1/1
feat1
#> An instance of class QFeatures (type: bulk) with 2 sets:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 3 rows and 2 columns 

## Aggregate peptides into proteins
feat1 <- aggregateFeatures(feat1, "peptides", "Protein", name = "proteins")
#> Aggregated: 1/1
feat1
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 3 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 2 rows and 2 columns 

assay(feat1[[1]])
#>       S1 S2
#> PSM1   1 11
#> PSM2   2 12
#> PSM3   3 13
#> PSM4   4 14
#> PSM5   5 15
#> PSM6   6 16
#> PSM7   7 17
#> PSM8   8 18
#> PSM9   9 19
#> PSM10 10 20
assay(feat1[[2]])
#>              S1   S2
#> ELGNDAYK    5.0 15.0
#> IAEESNFPFIK 8.5 18.5
#> SYGFNAAR    2.0 12.0
aggcounts(feat1[[2]])
#>             S1 S2
#> ELGNDAYK     3  3
#> IAEESNFPFIK  4  4
#> SYGFNAAR     3  3
assay(feat1[[3]])
#>        S1   S2
#> ProtA 3.5 13.5
#> ProtB 8.5 18.5
aggcounts(feat1[[3]])
#>       S1 S2
#> ProtA  2  2
#> ProtB  1  1

## --------------------------------------------
## Aggregation with missing quantitative values
## --------------------------------------------
data(ft_na)
ft_na
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] na: SummarizedExperiment with 4 rows and 3 columns 

assay(ft_na[[1]])
#>    A  B  C
#> a NA  5  9
#> b  2  6 10
#> c  3 NA 11
#> d NA  8 12
rowData(ft_na[[1]])
#> DataFrame with 4 rows and 2 columns
#>           X           Y
#>   <integer> <character>
#> a         1           A
#> b         2           B
#> c         1           A
#> d         2           B

## By default, missing values are propagated
ft2 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Aggregated: 1/1
assay(ft2[[2]])
#>    A  B  C
#> 1 NA NA 20
#> 2 NA 14 22
aggcounts(ft2[[2]])
#>   A B C
#> 1 1 1 2
#> 2 1 2 2

## The rowData .n variable tallies number of initial rows that
## were aggregated (irrespective of NAs) for all the samples.
rowData(ft2[[2]])
#> DataFrame with 2 rows and 3 columns
#>           X           Y        .n
#>   <integer> <character> <integer>
#> 1         1           A         2
#> 2         2           B         2

## Ignored when setting na.rm = TRUE
ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums, na.rm = TRUE)
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Aggregated: 1/1
assay(ft3[[2]])
#>   A  B  C
#> 1 3  5 20
#> 2 2 14 22
aggcounts(ft3[[2]])
#>   A B C
#> 1 1 1 2
#> 2 1 2 2

## -----------------------------------------------
## Aggregation with missing values in the row data
## -----------------------------------------------
## Row data results without any NAs, which includes the
## Y variables
rowData(ft2[[2]])
#> DataFrame with 2 rows and 3 columns
#>           X           Y        .n
#>   <integer> <character> <integer>
#> 1         1           A         2
#> 2         2           B         2

## Missing value in the Y feature variable
rowData(ft_na[[1]])[1, "Y"] <- NA
rowData(ft_na[[1]])
#> DataFrame with 4 rows and 2 columns
#>           X           Y
#>   <integer> <character>
#> a         1          NA
#> b         2           B
#> c         1           A
#> d         2           B

ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
#> Your quantitative and row data contain missing values. Please read the
#> relevant section(s) in the aggregateFeatures manual page regarding the
#> effects of missing values on data aggregation.
#> Aggregated: 1/1
## The Y feature variable has been dropped!
assay(ft3[[2]])
#>    A  B  C
#> 1 NA NA 20
#> 2 NA 14 22
rowData(ft3[[2]])
#> DataFrame with 2 rows and 2 columns
#>           X        .n
#>   <integer> <integer>
#> 1         1         2
#> 2         2         2

## --------------------------------------------
## Using a peptide-by-proteins adjacency matrix
## --------------------------------------------

## Let's use assay peptides from object feat1 and
## define that peptide SYGFNAAR maps to proteins
## Prot A and B

se <- feat1[["peptides"]]
rowData(se)$Protein[3] <- c("ProtA;ProtB")
rowData(se)
#> DataFrame with 3 rows and 4 columns
#>                  Sequence       Protein      location        .n
#>               <character>   <character>   <character> <integer>
#> ELGNDAYK         ELGNDAYK         ProtA Mitochondr...         3
#> IAEESNFPFIK IAEESNFPFI...         ProtB       unknown         4
#> SYGFNAAR         SYGFNAAR ProtA;Prot... Mitochondr...         3

## This can also be defined using anadjacency matrix, manual
## encoding here. See PSMatch::makeAdjacencyMatrix() for a
## function that does it automatically.
adj <- matrix(0, nrow = 3, ncol = 2,
              dimnames = list(rownames(se),
                              c("ProtA", "ProtB")))
adj[1, 1] <- adj[2, 2] <- adj[3, 1:2] <- 1
adj
#>             ProtA ProtB
#> ELGNDAYK        1     0
#> IAEESNFPFIK     0     1
#> SYGFNAAR        1     1

adjacencyMatrix(se) <- adj
rowData(se)
#> DataFrame with 3 rows and 5 columns
#>                  Sequence       Protein      location        .n adjacencyMatrix
#>               <character>   <character>   <character> <integer>     <dgCMatrix>
#> ELGNDAYK         ELGNDAYK         ProtA Mitochondr...         3             1:0
#> IAEESNFPFIK IAEESNFPFI...         ProtB       unknown         4             0:1
#> SYGFNAAR         SYGFNAAR ProtA;Prot... Mitochondr...         3             1:1
adjacencyMatrix(se)
#> 3 x 2 sparse Matrix of class "dgCMatrix"
#>             ProtA ProtB
#> ELGNDAYK        1     .
#> IAEESNFPFIK     .     1
#> SYGFNAAR        1     1

## Aggregation using the adjacency matrix
se2 <- aggregateFeatures(se, fcol = "adjacencyMatrix",
                         fun = MsCoreUtils::colMeansMat)

## Peptide SYGFNAAR was taken into account in both ProtA and ProtB
## aggregations.
assay(se2)
#>         S1    S2
#> ProtA 3.50 13.50
#> ProtB 5.25 15.25


## Aggregation by matrix on a QFeature object works as with a
## vector
ft <- QFeatures(list(peps = se))
ft <- aggregateFeatures(ft, "peps", "adjacencyMatrix", name = "protsByMat",
                        fun = MsCoreUtils::colMeansMat)
#> Aggregated: 1/1
assay(ft[[2]])
#>         S1    S2
#> ProtA 3.50 13.50
#> ProtB 5.25 15.25
rowData(ft[[2]])
#> DataFrame with 2 rows and 1 column
#>              .n
#>       <integer>
#> ProtA         2
#> ProtB         2
```
