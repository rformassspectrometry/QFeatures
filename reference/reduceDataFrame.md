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
#> 1           62      1.11         1           C      1-10     1        63
#> 2           25     -0.74         2           W      2-11     2        26
#> 3           96     -0.25         3           U      3-12     2        97
#> 4           91     -0.10         4           X      4-13     2        92
#> 5           40     -0.57         5           T      5-14     5        41
#> ...        ...       ...       ...         ...       ...   ...       ...
#> 996         22     -0.65       996           E  996-1005     1        23
#> 997         99     -0.10       997           V  997-1006     5       100
#> 998         63     -0.94       998           E  998-1007     4        64
#> 999        100      0.43       999           X  999-1008     4       101
#> 1000        25      0.27      1000           J 1000-1009     1        26

## Shinks the DataFrame
df2 <- reduceDataFrame(df, df$k)
df2
#> DataFrame with 100 rows and 7 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1 -0.17, 0.37,-1.05,...  58,193,323,...       C,F,G,...
#> 2           2  1.55,-0.08, 0.12,... 118,295,439,...       E,P,I,...
#> 3           3  0.29,-2.80, 1.29,...  32,140,477,...       X,E,U,...
#> 4           4    0.87,0.24,1.20,...  18,194,319,...       Q,Q,V,...
#> 5           5  0.86,-0.89,-0.28,...   15,63,121,...       J,Q,V,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.25,-0.17,-0.13,...   3,261,410,...       U,M,C,...
#> 97         97  0.31,-0.27,-0.57,...   37,42,104,...       L,C,C,...
#> 98         98  1.14, 0.84,-0.17,... 144,163,221,...       W,T,P,...
#> 99         99    0.04,2.16,0.47,...   33,82,159,...       L,L,A,...
#> 100       100  0.38, 0.27,-0.60,...   43,66,489,...       X,F,L,...
#>                              ir         r     invar
#>                   <IRangesList> <RleList> <numeric>
#> 1     58-67,193-202,323-332,... 2,4,2,...         2
#> 2   118-127,295-304,439-448,... 3,5,1,...         3
#> 3     32-41,140-149,477-486,... 4,2,4,...         4
#> 4     18-27,194-203,319-328,... 1,2,3,...         5
#> 5       15-24,63-72,121-130,... 2,5,3,...         6
#> ...                         ...       ...       ...
#> 96     3-12,261-270,410-419,... 2,4,1,...        97
#> 97      37-46,42-51,104-113,... 3,4,1,...        98
#> 98  144-153,163-172,221-230,... 4,1,4,...        99
#> 99      33-42,82-91,159-168,... 4,1,4,...       100
#> 100     43-52,66-75,489-498,... 1,1,4,...       101

## With a tally of the number of members in each group
reduceDataFrame(df, df$k, count = TRUE)
#> DataFrame with 100 rows and 8 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1 -0.17, 0.37,-1.05,...  58,193,323,...       C,F,G,...
#> 2           2  1.55,-0.08, 0.12,... 118,295,439,...       E,P,I,...
#> 3           3  0.29,-2.80, 1.29,...  32,140,477,...       X,E,U,...
#> 4           4    0.87,0.24,1.20,...  18,194,319,...       Q,Q,V,...
#> 5           5  0.86,-0.89,-0.28,...   15,63,121,...       J,Q,V,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.25,-0.17,-0.13,...   3,261,410,...       U,M,C,...
#> 97         97  0.31,-0.27,-0.57,...   37,42,104,...       L,C,C,...
#> 98         98  1.14, 0.84,-0.17,... 144,163,221,...       W,T,P,...
#> 99         99    0.04,2.16,0.47,...   33,82,159,...       L,L,A,...
#> 100       100  0.38, 0.27,-0.60,...   43,66,489,...       X,F,L,...
#>                              ir         r     invar        .n
#>                   <IRangesList> <RleList> <numeric> <integer>
#> 1     58-67,193-202,323-332,... 2,4,2,...         2        12
#> 2   118-127,295-304,439-448,... 3,5,1,...         3         9
#> 3     32-41,140-149,477-486,... 4,2,4,...         4         9
#> 4     18-27,194-203,319-328,... 1,2,3,...         5        12
#> 5       15-24,63-72,121-130,... 2,5,3,...         6         9
#> ...                         ...       ...       ...       ...
#> 96     3-12,261-270,410-419,... 2,4,1,...        97         7
#> 97      37-46,42-51,104-113,... 3,4,1,...        98        12
#> 98  144-153,163-172,221-230,... 4,1,4,...        99        10
#> 99      33-42,82-91,159-168,... 4,1,4,...       100        16
#> 100     43-52,66-75,489-498,... 1,1,4,...       101        11

## Much faster, but more crowded result
df3 <- reduceDataFrame(df, df$k, simplify = FALSE)
df3
#> DataFrame with 100 rows and 7 columns
#>                   k                     x               y               z
#>       <IntegerList>         <NumericList>   <IntegerList> <CharacterList>
#> 1         1,1,1,... -0.17, 0.37,-1.05,...  58,193,323,...       C,F,G,...
#> 2         2,2,2,...  1.55,-0.08, 0.12,... 118,295,439,...       E,P,I,...
#> 3         3,3,3,...  0.29,-2.80, 1.29,...  32,140,477,...       X,E,U,...
#> 4         4,4,4,...    0.87,0.24,1.20,...  18,194,319,...       Q,Q,V,...
#> 5         5,5,5,...  0.86,-0.89,-0.28,...   15,63,121,...       J,Q,V,...
#> ...             ...                   ...             ...             ...
#> 96     96,96,96,... -0.25,-0.17,-0.13,...   3,261,410,...       U,M,C,...
#> 97     97,97,97,...  0.31,-0.27,-0.57,...   37,42,104,...       L,C,C,...
#> 98     98,98,98,...  1.14, 0.84,-0.17,... 144,163,221,...       W,T,P,...
#> 99     99,99,99,...    0.04,2.16,0.47,...   33,82,159,...       L,L,A,...
#> 100 100,100,100,...  0.38, 0.27,-0.60,...   43,66,489,...       X,F,L,...
#>                              ir         r           invar
#>                   <IRangesList> <RleList>   <NumericList>
#> 1     58-67,193-202,323-332,... 2,4,2,...       2,2,2,...
#> 2   118-127,295-304,439-448,... 3,5,1,...       3,3,3,...
#> 3     32-41,140-149,477-486,... 4,2,4,...       4,4,4,...
#> 4     18-27,194-203,319-328,... 1,2,3,...       5,5,5,...
#> 5       15-24,63-72,121-130,... 2,5,3,...       6,6,6,...
#> ...                         ...       ...             ...
#> 96     3-12,261-270,410-419,... 2,4,1,...    97,97,97,...
#> 97      37-46,42-51,104-113,... 3,4,1,...    98,98,98,...
#> 98  144-153,163-172,221-230,... 4,1,4,...    99,99,99,...
#> 99      33-42,82-91,159-168,... 4,1,4,... 100,100,100,...
#> 100     43-52,66-75,489-498,... 1,1,4,... 101,101,101,...

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
