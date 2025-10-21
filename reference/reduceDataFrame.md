# Reduces and expands a `DataFrame`

A long dataframe can be *reduced* by mergeing certain rows into a single
one. These new variables are constructed as a `SimpleList` containing
all the original values. Invariant columns, i.e columns that have the
same value along all the rows that need to be merged, can be shrunk into
a new variables containing that invariant value (rather than in list
columns). The grouping of rows, i.e. the rows that need to be shrunk
together as one, is defined by a vector.

The opposite operation is *expand*. But note that for a `DataFrame` to
be expanded back, it must not to be simplified.

## Usage

``` r
reduceDataFrame(x, k, count = FALSE, simplify = TRUE, drop = FALSE)

expandDataFrame(x, k = NULL)
```

## Arguments

- x:

  The `DataFrame` to be reduced or expanded.

- k:

  A ‘vector’ of length `nrow(x)` defining the grouping based on which
  the `DataFrame` will be shrunk.

- count:

  `logical(1)` specifying of an additional column (called by default
  `.n`) with the tally of rows shrunk into on new row should be added.
  Note that if already existing, `.n` will be silently overwritten.

- simplify:

  A `logical(1)` defining if invariant columns should be converted to
  simple lists. Default is `TRUE`.

- drop:

  A `logical(1)` specifying whether the non-invariant columns should be
  dropped altogether. Default is `FALSE`.

## Value

An expanded (reduced) `DataFrame`.

## Missing values

Missing values do have an important effect on `reduce`. Unless all
values to be reduces are missing, they will result in an non-invariant
column, and will be dropped with `drop = TRUE`. See the example below.

The presence of missing values can have side effects in higher level
functions that rely on reduction of `DataFrame` objects.

## Author

Laurent Gatto

## Examples

``` r
library("IRanges")

k <- sample(100, 1e3, replace = TRUE)
df <- DataFrame(k = k,
                x = round(rnorm(length(k)), 2),
                y = seq_len(length(k)),
                z = sample(LETTERS, length(k), replace = TRUE),
                ir = IRanges(seq_along(k), width = 10),
                r = Rle(sample(5, length(k), replace = TRUE)),
                invar = k + 1)
df
#> DataFrame with 1000 rows and 7 columns
#>              k         x         y           z        ir     r     invar
#>      <integer> <numeric> <integer> <character> <IRanges> <Rle> <numeric>
#> 1            7      0.90         1           U      1-10     1         8
#> 2           39      1.39         2           I      2-11     4        40
#> 3           63     -0.82         3           M      3-12     5        64
#> 4           37     -0.58         4           O      4-13     1        38
#> 5           64     -1.23         5           K      5-14     5        65
#> ...        ...       ...       ...         ...       ...   ...       ...
#> 996         18      1.41       996           R  996-1005     5        19
#> 997         11      0.47       997           A  997-1006     3        12
#> 998         51      0.72       998           S  998-1007     5        52
#> 999         30     -1.26       999           X  999-1008     3        31
#> 1000        75     -1.38      1000           P 1000-1009     4        76

## Shinks the DataFrame
df2 <- reduceDataFrame(df, df$k)
df2
#> DataFrame with 100 rows and 7 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1 -2.16, 0.26,-0.04,...  33,168,298,...       E,E,T,...
#> 2           2  0.47,-0.23,-2.66,...  93,270,414,...       E,E,T,...
#> 3           3  2.09,-0.30, 0.85,...   7,115,452,...       P,P,F,...
#> 4           4  0.66,-0.75, 0.10,... 169,294,296,...       O,E,O,...
#> 5           5  1.21,-1.86, 1.83,...   38,96,109,...       C,X,W,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.84,-1.08, 1.44,... 236,385,388,...       T,G,A,...
#> 97         97  0.36,-1.79,-2.64,...    12,17,79,...       X,L,Q,...
#> 98         98 -1.08, 0.31, 0.84,... 119,138,196,...       I,R,L,...
#> 99         99  0.30,-0.22, 1.70,...    8,57,134,...       B,H,W,...
#> 100       100 -1.30,-0.24,-0.21,...   18,41,464,...       I,J,B,...
#>                              ir         r     invar
#>                   <IRangesList> <RleList> <numeric>
#> 1     33-42,168-177,298-307,... 2,3,5,...         2
#> 2    93-102,270-279,414-423,... 5,5,2,...         3
#> 3      7-16,115-124,452-461,... 2,2,1,...         4
#> 4   169-178,294-303,296-305,... 5,3,4,...         5
#> 5      38-47,96-105,109-118,... 1,5,5,...         6
#> ...                         ...       ...       ...
#> 96  236-245,385-394,388-397,... 2,5,3,...        97
#> 97        12-21,17-26,79-88,... 5,5,4,...        98
#> 98  119-128,138-147,196-205,... 2,4,5,...        99
#> 99       8-17,57-66,134-143,... 2,4,3,...       100
#> 100     18-27,41-50,464-473,... 1,2,5,...       101

## With a tally of the number of members in each group
reduceDataFrame(df, df$k, count = TRUE)
#> DataFrame with 100 rows and 8 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1 -2.16, 0.26,-0.04,...  33,168,298,...       E,E,T,...
#> 2           2  0.47,-0.23,-2.66,...  93,270,414,...       E,E,T,...
#> 3           3  2.09,-0.30, 0.85,...   7,115,452,...       P,P,F,...
#> 4           4  0.66,-0.75, 0.10,... 169,294,296,...       O,E,O,...
#> 5           5  1.21,-1.86, 1.83,...   38,96,109,...       C,X,W,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.84,-1.08, 1.44,... 236,385,388,...       T,G,A,...
#> 97         97  0.36,-1.79,-2.64,...    12,17,79,...       X,L,Q,...
#> 98         98 -1.08, 0.31, 0.84,... 119,138,196,...       I,R,L,...
#> 99         99  0.30,-0.22, 1.70,...    8,57,134,...       B,H,W,...
#> 100       100 -1.30,-0.24,-0.21,...   18,41,464,...       I,J,B,...
#>                              ir         r     invar        .n
#>                   <IRangesList> <RleList> <numeric> <integer>
#> 1     33-42,168-177,298-307,... 2,3,5,...         2        12
#> 2    93-102,270-279,414-423,... 5,5,2,...         3         9
#> 3      7-16,115-124,452-461,... 2,2,1,...         4         9
#> 4   169-178,294-303,296-305,... 5,3,4,...         5        11
#> 5      38-47,96-105,109-118,... 1,5,5,...         6         8
#> ...                         ...       ...       ...       ...
#> 96  236-245,385-394,388-397,... 2,5,3,...        97         6
#> 97        12-21,17-26,79-88,... 5,5,4,...        98        12
#> 98  119-128,138-147,196-205,... 2,4,5,...        99        10
#> 99       8-17,57-66,134-143,... 2,4,3,...       100        16
#> 100     18-27,41-50,464-473,... 1,2,5,...       101        11

## Much faster, but more crowded result
df3 <- reduceDataFrame(df, df$k, simplify = FALSE)
df3
#> DataFrame with 100 rows and 7 columns
#>                   k                     x               y               z
#>       <IntegerList>         <NumericList>   <IntegerList> <CharacterList>
#> 1         1,1,1,... -2.16, 0.26,-0.04,...  33,168,298,...       E,E,T,...
#> 2         2,2,2,...  0.47,-0.23,-2.66,...  93,270,414,...       E,E,T,...
#> 3         3,3,3,...  2.09,-0.30, 0.85,...   7,115,452,...       P,P,F,...
#> 4         4,4,4,...  0.66,-0.75, 0.10,... 169,294,296,...       O,E,O,...
#> 5         5,5,5,...  1.21,-1.86, 1.83,...   38,96,109,...       C,X,W,...
#> ...             ...                   ...             ...             ...
#> 96     96,96,96,... -0.84,-1.08, 1.44,... 236,385,388,...       T,G,A,...
#> 97     97,97,97,...  0.36,-1.79,-2.64,...    12,17,79,...       X,L,Q,...
#> 98     98,98,98,... -1.08, 0.31, 0.84,... 119,138,196,...       I,R,L,...
#> 99     99,99,99,...  0.30,-0.22, 1.70,...    8,57,134,...       B,H,W,...
#> 100 100,100,100,... -1.30,-0.24,-0.21,...   18,41,464,...       I,J,B,...
#>                              ir         r           invar
#>                   <IRangesList> <RleList>   <NumericList>
#> 1     33-42,168-177,298-307,... 2,3,5,...       2,2,2,...
#> 2    93-102,270-279,414-423,... 5,5,2,...       3,3,3,...
#> 3      7-16,115-124,452-461,... 2,2,1,...       4,4,4,...
#> 4   169-178,294-303,296-305,... 5,3,4,...       5,5,5,...
#> 5      38-47,96-105,109-118,... 1,5,5,...       6,6,6,...
#> ...                         ...       ...             ...
#> 96  236-245,385-394,388-397,... 2,5,3,...    97,97,97,...
#> 97        12-21,17-26,79-88,... 5,5,4,...    98,98,98,...
#> 98  119-128,138-147,196-205,... 2,4,5,...    99,99,99,...
#> 99       8-17,57-66,134-143,... 2,4,3,... 100,100,100,...
#> 100     18-27,41-50,464-473,... 1,2,5,... 101,101,101,...

## Drop all non-invariant columns
reduceDataFrame(df, df$k, drop = TRUE)
#> DataFrame with 100 rows and 2 columns
#>             k     invar
#>     <integer> <numeric>
#> 1           1         2
#> 2           2         3
#> 3           3         4
#> 4           4         5
#> 5           5         6
#> ...       ...       ...
#> 96         96        97
#> 97         97        98
#> 98         98        99
#> 99         99       100
#> 100       100       101

## Missing values
d <- DataFrame(k = rep(1:3, each = 3),
               x = letters[1:9],
               y = rep(letters[1:3], each = 3),
               y2 = rep(letters[1:3], each = 3))
d
#> DataFrame with 9 rows and 4 columns
#>           k           x           y          y2
#>   <integer> <character> <character> <character>
#> 1         1           a           a           a
#> 2         1           b           a           a
#> 3         1           c           a           a
#> 4         2           d           b           b
#> 5         2           e           b           b
#> 6         2           f           b           b
#> 7         3           g           c           c
#> 8         3           h           c           c
#> 9         3           i           c           c

## y is invariant and can be simplified
reduceDataFrame(d, d$k)
#> DataFrame with 3 rows and 4 columns
#>           k               x           y          y2
#>   <integer> <CharacterList> <character> <character>
#> 1         1           a,b,c           a           a
#> 2         2           d,e,f           b           b
#> 3         3           g,h,i           c           c
## y isn't not dropped
reduceDataFrame(d, d$k, drop = TRUE)
#> DataFrame with 3 rows and 3 columns
#>           k           y          y2
#>   <integer> <character> <character>
#> 1         1           a           a
#> 2         2           b           b
#> 3         3           c           c

## BUT with a missing value
d[1, "y"] <- NA
d
#> DataFrame with 9 rows and 4 columns
#>           k           x           y          y2
#>   <integer> <character> <character> <character>
#> 1         1           a          NA           a
#> 2         1           b           a           a
#> 3         1           c           a           a
#> 4         2           d           b           b
#> 5         2           e           b           b
#> 6         2           f           b           b
#> 7         3           g           c           c
#> 8         3           h           c           c
#> 9         3           i           c           c

## y isn't invariant/simplified anymore
reduceDataFrame(d, d$k)
#> DataFrame with 3 rows and 4 columns
#>           k               x               y          y2
#>   <integer> <CharacterList> <CharacterList> <character>
#> 1         1           a,b,c          NA,a,a           a
#> 2         2           d,e,f           b,b,b           b
#> 3         3           g,h,i           c,c,c           c
## y now gets dropped
reduceDataFrame(d, d$k, drop = TRUE)
#> DataFrame with 3 rows and 2 columns
#>           k          y2
#>   <integer> <character>
#> 1         1           a
#> 2         2           b
#> 3         3           c
```
