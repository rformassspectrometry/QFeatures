# Filter features based on their rowData

The `filterFeatures` methods enables users to filter features based on a
variable in their `rowData`. The features matching the filter will be
returned as a new object of class `QFeatures`. The filters can be
provided as instances of class `AnnotationFilter` (see below) or as
formulas.

## Usage

``` r
VariableFilter(field, value, condition = "==", not = FALSE)

# S4 method for class 'QFeatures,AnnotationFilter'
filterFeatures(object, filter, i, na.rm = FALSE, keep = FALSE, ...)

# S4 method for class 'QFeatures,formula'
filterFeatures(object, filter, i, na.rm = FALSE, keep = FALSE, ...)

isDuplicated(x)
```

## Arguments

- field:

  `character(1)` refering to the name of the variable to apply the
  filter on.

- value:

  [`character()`](https://rdrr.io/r/base/character.html) or
  [`integer()`](https://rdrr.io/r/base/integer.html) value for the
  `CharacterVariableFilter` and `NumericVariableFilter` filters
  respectively.

- condition:

  `character(1)` defining the condition to be used in the filter. For
  `NumericVariableFilter`, one of `"=="`, `"!="`, `">"`, `"<"`, `">="`
  or `"<="`. For `CharacterVariableFilter`, one of `"=="`, `"!="`,
  `"startsWith"`, `"endsWith"` or `"contains"`. Default condition is
  `"=="`.

- not:

  `logical(1)` indicating whether the filtering should be negated or
  not. `TRUE` indicates is negated (!). `FALSE` indicates not negated.
  Default `not` is `FALSE`, so no negation.

- object:

  An instance of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md).

- filter:

  Either an instance of class AnnotationFilter or a formula.

- i:

  A numeric, logical or character vector pointing to the assay(s) to be
  filtered.

- na.rm:

  `logical(1)` indicating whether missing values should be removed.
  Default is `FALSE`.

- keep:

  `logical(1)` indicating whether to keep the features of assays for
  which at least one of the filtering variables are missing in the
  rowData. When `FALSE` (default), all such assay will contain 0
  features; when `TRUE`, the assays are untouched.

- ...:

  Additional parameters. Currently ignored.

- x:

  A [`vector()`](https://rdrr.io/r/base/vector.html) that will be
  checked for duplications.

## Value

An filtered `QFeature` object.

## The filtering procedure

`filterFeatures()` will go through each assay of the `QFeatures` object
and apply the filtering on the corresponding `rowData`. Features that do
not pass the filter condition are removed from the assay. In some cases,
one may want to filter for a variable present in some assay, but not in
other. There are two options: either provide `keep = FALSE` to remove
all features for those assays (and thus leaving an empty assay), or
provide `keep = TRUE` to ignore filtering for those assays.

Because features in a `QFeatures` object are linked between different
assays with `AssayLinks`, the links are automatically updated. However,
note that the function doesn't propagate the filter to parent assays.
For example, suppose a peptide assay with 4 peptides is linked to a
protein assay with 2 proteins (2 peptides mapped per protein) and you
apply `filterFeatures()`. All features pass the filter except for one
protein. The peptides mapped to that protein will remain in the
`QFeatures` object. If propagation of the filtering rules to parent
assay is desired, you may want to use `x[i, , ]` instead (see the
*Subsetting* section in `?QFeature`).

## Variable filters

The variable filters are filters as defined in the AnnotationFilter
package. In addition to the pre-defined filter, users can arbitrarily
set a field on which to operate. These arbitrary filters operate either
on a character variables (as `CharacterVariableFilter` objects) or
numerics (as `NumericVariableFilters` objects), which can be created
with the `VariableFilter` constructor.

## Helper functions

- The `isDuplicated()` function takes a vector (or rowData variable when
  used to filter features) as input, and return a logical of the same
  length, with elements set to `TRUE` for unique occurence, and `FALSE`
  otherwise. This function is different from
  [`duplicated()`](https://rdrr.io/r/base/duplicated.html), as here even
  the first occurence is set to `FALSE`. See
  [`createPrecursorId()`](https://rformassspectrometry.github.io/QFeatures/reference/createPrecursorId.md)
  for an application.

## See also

The
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
man page for subsetting and the `QFeatures` vignette provides an
extended example.

## Author

Laurent Gatto

## Examples

``` r

## ----------------------------------------
## Creating character and numberic
## variable filters
## ----------------------------------------

VariableFilter(field = "my_var",
               value = "value_to_keep",
               condition = "==")
#> class: CharacterVariableFilter 
#> condition: == 
#> value: value_to_keep 

VariableFilter(field = "my_num_var",
               value = 0.05,
               condition = "<=")
#> class: NumericVariableFilter 
#> condition: <= 
#> value: 0.05 

example(aggregateFeatures)
#> 
#> aggrgF> ## ---------------------------------------
#> aggrgF> ## An example QFeatures with PSM-level data
#> aggrgF> ## ---------------------------------------
#> aggrgF> data(feat1)
#> 
#> aggrgF> feat1
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 
#> 
#> aggrgF> ## Aggregate PSMs into peptides
#> aggrgF> feat1 <- aggregateFeatures(feat1, "psms", "Sequence", name = "peptides")
#> Aggregated: 1/1
#> 
#> aggrgF> feat1
#> An instance of class QFeatures (type: bulk) with 2 sets:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 3 rows and 2 columns 
#> 
#> aggrgF> ## Aggregate peptides into proteins
#> aggrgF> feat1 <- aggregateFeatures(feat1, "peptides", "Protein", name = "proteins")
#> Aggregated: 1/1
#> 
#> aggrgF> feat1
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 3 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 2 rows and 2 columns 
#> 
#> aggrgF> assay(feat1[[1]])
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
#> 
#> aggrgF> assay(feat1[[2]])
#>              S1   S2
#> ELGNDAYK    5.0 15.0
#> IAEESNFPFIK 8.5 18.5
#> SYGFNAAR    2.0 12.0
#> 
#> aggrgF> aggcounts(feat1[[2]])
#>             S1 S2
#> ELGNDAYK     3  3
#> IAEESNFPFIK  4  4
#> SYGFNAAR     3  3
#> 
#> aggrgF> assay(feat1[[3]])
#>        S1   S2
#> ProtA 3.5 13.5
#> ProtB 8.5 18.5
#> 
#> aggrgF> aggcounts(feat1[[3]])
#>       S1 S2
#> ProtA  2  2
#> ProtB  1  1
#> 
#> aggrgF> ## --------------------------------------------
#> aggrgF> ## Aggregation with missing quantitative values
#> aggrgF> ## --------------------------------------------
#> aggrgF> data(ft_na)
#> 
#> aggrgF> ft_na
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] na: SummarizedExperiment with 4 rows and 3 columns 
#> 
#> aggrgF> assay(ft_na[[1]])
#>    A  B  C
#> a NA  5  9
#> b  2  6 10
#> c  3 NA 11
#> d NA  8 12
#> 
#> aggrgF> rowData(ft_na[[1]])
#> DataFrame with 4 rows and 2 columns
#>           X           Y
#>   <integer> <character>
#> a         1           A
#> b         2           B
#> c         1           A
#> d         2           B
#> 
#> aggrgF> ## By default, missing values are propagated
#> aggrgF> ft2 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Aggregated: 1/1
#> 
#> aggrgF> assay(ft2[[2]])
#>    A  B  C
#> 1 NA NA 20
#> 2 NA 14 22
#> 
#> aggrgF> aggcounts(ft2[[2]])
#>   A B C
#> 1 1 1 2
#> 2 1 2 2
#> 
#> aggrgF> ## The rowData .n variable tallies number of initial rows that
#> aggrgF> ## were aggregated (irrespective of NAs) for all the samples.
#> aggrgF> rowData(ft2[[2]])
#> DataFrame with 2 rows and 3 columns
#>           X           Y        .n
#>   <integer> <character> <integer>
#> 1         1           A         2
#> 2         2           B         2
#> 
#> aggrgF> ## Ignored when setting na.rm = TRUE
#> aggrgF> ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums, na.rm = TRUE)
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Aggregated: 1/1
#> 
#> aggrgF> assay(ft3[[2]])
#>   A  B  C
#> 1 3  5 20
#> 2 2 14 22
#> 
#> aggrgF> aggcounts(ft3[[2]])
#>   A B C
#> 1 1 1 2
#> 2 1 2 2
#> 
#> aggrgF> ## -----------------------------------------------
#> aggrgF> ## Aggregation with missing values in the row data
#> aggrgF> ## -----------------------------------------------
#> aggrgF> ## Row data results without any NAs, which includes the
#> aggrgF> ## Y variables
#> aggrgF> rowData(ft2[[2]])
#> DataFrame with 2 rows and 3 columns
#>           X           Y        .n
#>   <integer> <character> <integer>
#> 1         1           A         2
#> 2         2           B         2
#> 
#> aggrgF> ## Missing value in the Y feature variable
#> aggrgF> rowData(ft_na[[1]])[1, "Y"] <- NA
#> 
#> aggrgF> rowData(ft_na[[1]])
#> DataFrame with 4 rows and 2 columns
#>           X           Y
#>   <integer> <character>
#> a         1          NA
#> b         2           B
#> c         1           A
#> d         2           B
#> 
#> aggrgF> ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
#> Your quantitative and row data contain missing values. Please read the
#> relevant section(s) in the aggregateFeatures manual page regarding the
#> effects of missing values on data aggregation.
#> Aggregated: 1/1
#> 
#> aggrgF> ## The Y feature variable has been dropped!
#> aggrgF> assay(ft3[[2]])
#>    A  B  C
#> 1 NA NA 20
#> 2 NA 14 22
#> 
#> aggrgF> rowData(ft3[[2]])
#> DataFrame with 2 rows and 2 columns
#>           X        .n
#>   <integer> <integer>
#> 1         1         2
#> 2         2         2
#> 
#> aggrgF> ## --------------------------------------------
#> aggrgF> ## Using a peptide-by-proteins adjacency matrix
#> aggrgF> ## --------------------------------------------
#> aggrgF> 
#> aggrgF> ## Let's use assay peptides from object feat1 and
#> aggrgF> ## define that peptide SYGFNAAR maps to proteins
#> aggrgF> ## Prot A and B
#> aggrgF> 
#> aggrgF> se <- feat1[["peptides"]]
#> 
#> aggrgF> rowData(se)$Protein[3] <- c("ProtA;ProtB")
#> 
#> aggrgF> rowData(se)
#> DataFrame with 3 rows and 4 columns
#>                  Sequence       Protein      location        .n
#>               <character>   <character>   <character> <integer>
#> ELGNDAYK         ELGNDAYK         ProtA Mitochondr...         3
#> IAEESNFPFIK IAEESNFPFI...         ProtB       unknown         4
#> SYGFNAAR         SYGFNAAR ProtA;Prot... Mitochondr...         3
#> 
#> aggrgF> ## This can also be defined using anadjacency matrix, manual
#> aggrgF> ## encoding here. See PSMatch::makeAdjacencyMatrix() for a
#> aggrgF> ## function that does it automatically.
#> aggrgF> adj <- matrix(0, nrow = 3, ncol = 2,
#> aggrgF+               dimnames = list(rownames(se),
#> aggrgF+                               c("ProtA", "ProtB")))
#> 
#> aggrgF> adj[1, 1] <- adj[2, 2] <- adj[3, 1:2] <- 1
#> 
#> aggrgF> adj
#>             ProtA ProtB
#> ELGNDAYK        1     0
#> IAEESNFPFIK     0     1
#> SYGFNAAR        1     1
#> 
#> aggrgF> adjacencyMatrix(se) <- adj
#> 
#> aggrgF> rowData(se)
#> DataFrame with 3 rows and 5 columns
#>                  Sequence       Protein      location        .n adjacencyMatrix
#>               <character>   <character>   <character> <integer>     <dgCMatrix>
#> ELGNDAYK         ELGNDAYK         ProtA Mitochondr...         3             1:0
#> IAEESNFPFIK IAEESNFPFI...         ProtB       unknown         4             0:1
#> SYGFNAAR         SYGFNAAR ProtA;Prot... Mitochondr...         3             1:1
#> 
#> aggrgF> adjacencyMatrix(se)
#> 3 x 2 sparse Matrix of class "dgCMatrix"
#>             ProtA ProtB
#> ELGNDAYK        1     .
#> IAEESNFPFIK     .     1
#> SYGFNAAR        1     1
#> 
#> aggrgF> ## Aggregation using the adjacency matrix
#> aggrgF> se2 <- aggregateFeatures(se, fcol = "adjacencyMatrix",
#> aggrgF+                          fun = MsCoreUtils::colMeansMat)
#> 
#> aggrgF> ## Peptide SYGFNAAR was taken into account in both ProtA and ProtB
#> aggrgF> ## aggregations.
#> aggrgF> assay(se2)
#>         S1    S2
#> ProtA 3.50 13.50
#> ProtB 5.25 15.25
#> 
#> aggrgF> ## Aggregation by matrix on a QFeature object works as with a
#> aggrgF> ## vector
#> aggrgF> ft <- QFeatures(list(peps = se))
#> 
#> aggrgF> ft <- aggregateFeatures(ft, "peps", "adjacencyMatrix", name = "protsByMat",
#> aggrgF+                         fun = MsCoreUtils::colMeansMat)
#> Aggregated: 1/1
#> 
#> aggrgF> assay(ft[[2]])
#>         S1    S2
#> ProtA 3.50 13.50
#> ProtB 5.25 15.25
#> 
#> aggrgF> rowData(ft[[2]])
#> DataFrame with 2 rows and 1 column
#>              .n
#>       <integer>
#> ProtA         2
#> ProtB         2

## ----------------------------------------------------------------
## Filter all features that are associated to the Mitochondrion in
## the location feature variable. This variable is present in all
## assays.
## ----------------------------------------------------------------

## using the forumla interface, exact mathc
filterFeatures(feat1, ~  location == "Mitochondrion")
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 

## using the forumula intefrace, martial match
filterFeatures(feat1, ~startsWith(location, "Mito"))
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 

## using a user-defined character filter
filterFeatures(feat1, VariableFilter("location", "Mitochondrion"))
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 

## using a user-defined character filter with partial match
filterFeatures(feat1, VariableFilter("location", "Mito", "startsWith"))
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 
filterFeatures(feat1, VariableFilter("location", "itochon", "contains"))
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 

## ----------------------------------------------------------------
## Filter all features that aren't marked as unknown (sub-cellular
## location) in the feature variable
## ----------------------------------------------------------------

## using a user-defined character filter
filterFeatures(feat1, VariableFilter("location", "unknown", condition = "!="))
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 

## using the forumula interface
filterFeatures(feat1, ~ location != "unknown")
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 

## ----------------------------------------------------------------
## Filter features that have a p-values lower or equal to 0.03
## ----------------------------------------------------------------

## using a user-defined numeric filter
filterFeatures(feat1, VariableFilter("pval", 0.03, "<="))
#> 'pval' found in 1 out of 3 assay(s).
#> No filter applied to the following assay(s) because one or more
#> filtering variables are missing in the rowData: peptides, proteins. You
#> can control whether to remove or keep the features using the 'keep'
#> argument (see '?filterFeature').
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 3 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 0 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 0 rows and 2 columns 

## using the formula interface
filterFeatures(feat1, ~ pval <= 0.03)
#> 'pval' found in 1 out of 3 assay(s).
#> No filter applied to the following assay(s) because one or more
#> filtering variables are missing in the rowData: peptides, proteins. You
#> can control whether to remove or keep the features using the 'keep'
#> argument (see '?filterFeature').
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 3 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 0 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 0 rows and 2 columns 

## you can also remove all p-values that are NA (if any)
filterFeatures(feat1, ~ !is.na(pval))
#> 'pval' found in 1 out of 3 assay(s).
#> No filter applied to the following assay(s) because one or more
#> filtering variables are missing in the rowData: peptides, proteins. You
#> can control whether to remove or keep the features using the 'keep'
#> argument (see '?filterFeature').
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 0 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 0 rows and 2 columns 

## ----------------------------------------------------------------
## Negative control - filtering for an non-existing markers value,
## returning empty results.
## ----------------------------------------------------------------

filterFeatures(feat1, VariableFilter("location", "not"))
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 0 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 0 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 0 rows and 2 columns 

filterFeatures(feat1, ~ location == "not")
#> 'location' found in 3 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 0 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 0 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 0 rows and 2 columns 

## ----------------------------------------------------------------
## Filtering for a  missing feature variable. The outcome is controled
## by keep
## ----------------------------------------------------------------
data(feat2)

filterFeatures(feat2, ~ y < 0)
#> 'y' found in 2 out of 3 assay(s).
#> No filter applied to the following assay(s) because one or more
#> filtering variables are missing in the rowData: assay1. You can control
#> whether to remove or keep the features using the 'keep' argument (see
#> '?filterFeature').
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] assay1: SummarizedExperiment with 0 rows and 4 columns 
#>  [2] assay2: SummarizedExperiment with 1 rows and 4 columns 
#>  [3] assay3: SummarizedExperiment with 5 rows and 4 columns 

filterFeatures(feat2, ~ y < 0, keep = TRUE)
#> 'y' found in 2 out of 3 assay(s).
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] assay1: SummarizedExperiment with 10 rows and 4 columns 
#>  [2] assay2: SummarizedExperiment with 1 rows and 4 columns 
#>  [3] assay3: SummarizedExperiment with 5 rows and 4 columns 

## ----------------------------------------------------------------
## Example with missing values
## ----------------------------------------------------------------

data(feat1)
rowData(feat1[[1]])[1, "location"] <- NA
rowData(feat1[[1]])
#> DataFrame with 10 rows and 5 columns
#>            Sequence     Protein       Var      location      pval
#>         <character> <character> <integer>   <character> <numeric>
#> PSM1       SYGFNAAR       ProtA         1            NA     0.084
#> PSM2       SYGFNAAR       ProtA         2 Mitochondr...     0.077
#> PSM3       SYGFNAAR       ProtA         3 Mitochondr...     0.063
#> PSM4       ELGNDAYK       ProtA         4 Mitochondr...     0.073
#> PSM5       ELGNDAYK       ProtA         5 Mitochondr...     0.012
#> PSM6       ELGNDAYK       ProtA         6 Mitochondr...     0.011
#> PSM7  IAEESNFPFI...       ProtB         7       unknown     0.075
#> PSM8  IAEESNFPFI...       ProtB         8       unknown     0.038
#> PSM9  IAEESNFPFI...       ProtB         9       unknown     0.028
#> PSM10 IAEESNFPFI...       ProtB        10       unknown     0.097

## The row with the NA is not removed
rowData(filterFeatures(feat1, ~ location == "Mitochondrion")[[1]])
#> 'location' found in 1 out of 1 assay(s).
#> DataFrame with 6 rows and 5 columns
#>         Sequence     Protein       Var      location      pval
#>      <character> <character> <integer>   <character> <numeric>
#> PSM1    SYGFNAAR       ProtA         1            NA     0.084
#> PSM2    SYGFNAAR       ProtA         2 Mitochondr...     0.077
#> PSM3    SYGFNAAR       ProtA         3 Mitochondr...     0.063
#> PSM4    ELGNDAYK       ProtA         4 Mitochondr...     0.073
#> PSM5    ELGNDAYK       ProtA         5 Mitochondr...     0.012
#> PSM6    ELGNDAYK       ProtA         6 Mitochondr...     0.011
rowData(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = FALSE)[[1]])
#> 'location' found in 1 out of 1 assay(s).
#> DataFrame with 6 rows and 5 columns
#>         Sequence     Protein       Var      location      pval
#>      <character> <character> <integer>   <character> <numeric>
#> PSM1    SYGFNAAR       ProtA         1            NA     0.084
#> PSM2    SYGFNAAR       ProtA         2 Mitochondr...     0.077
#> PSM3    SYGFNAAR       ProtA         3 Mitochondr...     0.063
#> PSM4    ELGNDAYK       ProtA         4 Mitochondr...     0.073
#> PSM5    ELGNDAYK       ProtA         5 Mitochondr...     0.012
#> PSM6    ELGNDAYK       ProtA         6 Mitochondr...     0.011

## The row with the NA is removed
rowData(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = TRUE)[[1]])
#> 'location' found in 1 out of 1 assay(s).
#> DataFrame with 5 rows and 5 columns
#>         Sequence     Protein       Var      location      pval
#>      <character> <character> <integer>   <character> <numeric>
#> PSM2    SYGFNAAR       ProtA         2 Mitochondr...     0.077
#> PSM3    SYGFNAAR       ProtA         3 Mitochondr...     0.063
#> PSM4    ELGNDAYK       ProtA         4 Mitochondr...     0.073
#> PSM5    ELGNDAYK       ProtA         5 Mitochondr...     0.012
#> PSM6    ELGNDAYK       ProtA         6 Mitochondr...     0.011

## Note that in situations with missing values, it is possible to
## use the `%in%` operator or filter missing values out
## explicitly.

rowData(filterFeatures(feat1, ~ location %in% "Mitochondrion")[[1]])
#> 'location' found in 1 out of 1 assay(s).
#> DataFrame with 5 rows and 5 columns
#>         Sequence     Protein       Var      location      pval
#>      <character> <character> <integer>   <character> <numeric>
#> PSM2    SYGFNAAR       ProtA         2 Mitochondr...     0.077
#> PSM3    SYGFNAAR       ProtA         3 Mitochondr...     0.063
#> PSM4    ELGNDAYK       ProtA         4 Mitochondr...     0.073
#> PSM5    ELGNDAYK       ProtA         5 Mitochondr...     0.012
#> PSM6    ELGNDAYK       ProtA         6 Mitochondr...     0.011
rowData(filterFeatures(feat1, ~ location %in% c(NA, "Mitochondrion"))[[1]])
#> 'location' found in 1 out of 1 assay(s).
#> DataFrame with 6 rows and 5 columns
#>         Sequence     Protein       Var      location      pval
#>      <character> <character> <integer>   <character> <numeric>
#> PSM1    SYGFNAAR       ProtA         1            NA     0.084
#> PSM2    SYGFNAAR       ProtA         2 Mitochondr...     0.077
#> PSM3    SYGFNAAR       ProtA         3 Mitochondr...     0.063
#> PSM4    ELGNDAYK       ProtA         4 Mitochondr...     0.073
#> PSM5    ELGNDAYK       ProtA         5 Mitochondr...     0.012
#> PSM6    ELGNDAYK       ProtA         6 Mitochondr...     0.011

## Explicit handling
filterFeatures(feat1, ~ !is.na(location) & location == "Mitochondrion")
#> 'location' found in 1 out of 1 assay(s).
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] psms: SummarizedExperiment with 5 rows and 2 columns 

## Using the pipe operator
feat1 |>
   filterFeatures( ~ !is.na(location)) |>
   filterFeatures( ~ location == "Mitochondrion")
#> 'location' found in 1 out of 1 assay(s).
#> 'location' found in 1 out of 1 assay(s).
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] psms: SummarizedExperiment with 5 rows and 2 columns 
```
