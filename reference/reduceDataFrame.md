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
#> 1           63      1.39         1           M      1-10     5        64
#> 2           37     -0.82         2           O      2-11     1        38
#> 3           64     -0.58         3           K      3-12     5        65
#> 4           18     -1.23         4           V      4-13     5        19
#> 5            3      1.28         5           P      5-14     2         4
#> ...        ...       ...       ...         ...       ...   ...       ...
#> 996         51      0.47       996           S  996-1005     5        52
#> 997         30      0.72       997           X  997-1006     3        31
#> 998         75     -1.26       998           P  998-1007     4        76
#> 999          1     -1.38       999           Q  999-1008     2         2
#> 1000        13     -0.53      1000           T 1000-1009     3        14

## Shinks the DataFrame
df2 <- reduceDataFrame(df, df$k)
df2
#> DataFrame with 100 rows and 7 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1  0.87,-0.01, 1.02,...  31,166,296,...       E,E,T,...
#> 2           2 -1.78, 0.00, 0.03,...  91,268,412,...       E,E,T,...
#> 3           3  1.28,-0.35, 0.47,...   5,113,450,...       P,P,F,...
#> 4           4  0.26,-0.32, 1.20,... 167,292,294,...       O,E,O,...
#> 5           5    1.10,0.38,0.00,...   36,94,107,...       C,X,W,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.87,-0.40, 1.24,... 234,383,386,...       T,G,A,...
#> 97         97 -0.90,-0.15, 0.40,...    10,15,77,...       X,L,Q,...
#> 98         98 -0.41, 1.28,-0.23,... 117,136,194,...       I,R,L,...
#> 99         99  2.09,-0.83,-0.26,...    6,55,132,...       B,H,W,...
#> 100       100 -1.79,-0.84, 1.91,...   16,39,462,...       I,J,B,...
#>                              ir         r     invar
#>                   <IRangesList> <RleList> <numeric>
#> 1     31-40,166-175,296-305,... 2,3,5,...         2
#> 2    91-100,268-277,412-421,... 5,5,2,...         3
#> 3      5-14,113-122,450-459,... 2,2,1,...         4
#> 4   167-176,292-301,294-303,... 5,3,4,...         5
#> 5      36-45,94-103,107-116,... 1,5,5,...         6
#> ...                         ...       ...       ...
#> 96  234-243,383-392,386-395,... 2,5,3,...        97
#> 97        10-19,15-24,77-86,... 5,5,4,...        98
#> 98  117-126,136-145,194-203,... 2,4,5,...        99
#> 99       6-15,55-64,132-141,... 2,4,3,...       100
#> 100     16-25,39-48,462-471,... 1,2,5,...       101

## With a tally of the number of members in each group
reduceDataFrame(df, df$k, count = TRUE)
#> DataFrame with 100 rows and 8 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1  0.87,-0.01, 1.02,...  31,166,296,...       E,E,T,...
#> 2           2 -1.78, 0.00, 0.03,...  91,268,412,...       E,E,T,...
#> 3           3  1.28,-0.35, 0.47,...   5,113,450,...       P,P,F,...
#> 4           4  0.26,-0.32, 1.20,... 167,292,294,...       O,E,O,...
#> 5           5    1.10,0.38,0.00,...   36,94,107,...       C,X,W,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.87,-0.40, 1.24,... 234,383,386,...       T,G,A,...
#> 97         97 -0.90,-0.15, 0.40,...    10,15,77,...       X,L,Q,...
#> 98         98 -0.41, 1.28,-0.23,... 117,136,194,...       I,R,L,...
#> 99         99  2.09,-0.83,-0.26,...    6,55,132,...       B,H,W,...
#> 100       100 -1.79,-0.84, 1.91,...   16,39,462,...       I,J,B,...
#>                              ir         r     invar        .n
#>                   <IRangesList> <RleList> <numeric> <integer>
#> 1     31-40,166-175,296-305,... 2,3,5,...         2        13
#> 2    91-100,268-277,412-421,... 5,5,2,...         3         9
#> 3      5-14,113-122,450-459,... 2,2,1,...         4         9
#> 4   167-176,292-301,294-303,... 5,3,4,...         5        11
#> 5      36-45,94-103,107-116,... 1,5,5,...         6         8
#> ...                         ...       ...       ...       ...
#> 96  234-243,383-392,386-395,... 2,5,3,...        97         6
#> 97        10-19,15-24,77-86,... 5,5,4,...        98        12
#> 98  117-126,136-145,194-203,... 2,4,5,...        99        10
#> 99       6-15,55-64,132-141,... 2,4,3,...       100        16
#> 100     16-25,39-48,462-471,... 1,2,5,...       101        11

## Much faster, but more crowded result
df3 <- reduceDataFrame(df, df$k, simplify = FALSE)
df3
#> DataFrame with 100 rows and 7 columns
#>                   k                     x               y               z
#>       <IntegerList>         <NumericList>   <IntegerList> <CharacterList>
#> 1         1,1,1,...  0.87,-0.01, 1.02,...  31,166,296,...       E,E,T,...
#> 2         2,2,2,... -1.78, 0.00, 0.03,...  91,268,412,...       E,E,T,...
#> 3         3,3,3,...  1.28,-0.35, 0.47,...   5,113,450,...       P,P,F,...
#> 4         4,4,4,...  0.26,-0.32, 1.20,... 167,292,294,...       O,E,O,...
#> 5         5,5,5,...    1.10,0.38,0.00,...   36,94,107,...       C,X,W,...
#> ...             ...                   ...             ...             ...
#> 96     96,96,96,... -0.87,-0.40, 1.24,... 234,383,386,...       T,G,A,...
#> 97     97,97,97,... -0.90,-0.15, 0.40,...    10,15,77,...       X,L,Q,...
#> 98     98,98,98,... -0.41, 1.28,-0.23,... 117,136,194,...       I,R,L,...
#> 99     99,99,99,...  2.09,-0.83,-0.26,...    6,55,132,...       B,H,W,...
#> 100 100,100,100,... -1.79,-0.84, 1.91,...   16,39,462,...       I,J,B,...
#>                              ir         r           invar
#>                   <IRangesList> <RleList>   <NumericList>
#> 1     31-40,166-175,296-305,... 2,3,5,...       2,2,2,...
#> 2    91-100,268-277,412-421,... 5,5,2,...       3,3,3,...
#> 3      5-14,113-122,450-459,... 2,2,1,...       4,4,4,...
#> 4   167-176,292-301,294-303,... 5,3,4,...       5,5,5,...
#> 5      36-45,94-103,107-116,... 1,5,5,...       6,6,6,...
#> ...                         ...       ...             ...
#> 96  234-243,383-392,386-395,... 2,5,3,...    97,97,97,...
#> 97        10-19,15-24,77-86,... 5,5,4,...    98,98,98,...
#> 98  117-126,136-145,194-203,... 2,4,5,...    99,99,99,...
#> 99       6-15,55-64,132-141,... 2,4,3,... 100,100,100,...
#> 100     16-25,39-48,462-471,... 1,2,5,... 101,101,101,...

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
