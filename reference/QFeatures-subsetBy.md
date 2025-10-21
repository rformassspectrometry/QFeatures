# Subset by feature name

This function will find the assays and features that match directly (by
name) or indirectly (through aggregation) the feature name.

The `subsetByFeature` function will first identify the assay that
contains the feature(s) `i` and filter the rows matching these feature
names exactly. It will then find, in the other assays, the features that
produces `i` through aggregation with the `aggregateQFeatures` function.

See
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
for an example.

## Arguments

- x:

  An instance of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md).

- y:

  A `character` of feature names present in an assay in `x`.

- ...:

  Additional parameters. Ignored.

## Value

An new instance of class
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
containing relevant assays and features.

## Examples

``` r
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

## Retrieve protein 'ProtA' and its 2 peptides and 6 PSMs
feat1["ProtA", , ]
#> An instance of class QFeatures (type: bulk) with 3 sets:
#> 
#>  [1] psms: SummarizedExperiment with 6 rows and 2 columns 
#>  [2] peptides: SummarizedExperiment with 2 rows and 2 columns 
#>  [3] proteins: SummarizedExperiment with 1 rows and 2 columns 
```
