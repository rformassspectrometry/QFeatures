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
#> 1           37     -0.31         1           M      1-10     5        38
#> 2           64      0.87         2           O      2-11     1        65
#> 3           18     -0.51         3           K      3-12     5        19
#> 4            3     -0.77         4           V      4-13     5         4
#> 5           99     -1.24         5           P      5-14     2       100
#> ...        ...       ...       ...         ...       ...   ...       ...
#> 996         30     -0.20       996           S  996-1005     5        31
#> 997         75     -0.66       997           X  997-1006     3        76
#> 998          1      0.45       998           P  998-1007     4         2
#> 999         13     -0.76       999           Q  999-1008     2        14
#> 1000        23     -0.14      1000           T 1000-1009     3        24

## Shinks the DataFrame
df2 <- reduceDataFrame(df, df$k)
df2
#> DataFrame with 100 rows and 7 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1 -0.69,-0.15,-1.90,...  30,165,295,...       P,I,M,...
#> 2           2  0.33, 1.16,-0.88,...  90,267,411,...       P,S,W,...
#> 3           3 -0.77,-0.14, 0.65,...   4,112,449,...       V,M,T,...
#> 4           4    0.73,0.52,1.67,... 166,291,293,...       E,U,H,...
#> 5           5 -0.47,-0.38,-0.10,...   35,93,106,...       K,C,P,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.32,-0.63, 1.09,... 233,382,385,...       U,U,B,...
#> 97         97  0.79, 0.43,-0.66,...     9,14,76,...       I,T,Z,...
#> 98         98  0.83,-0.53, 0.25,... 116,135,193,...       A,C,S,...
#> 99         99 -1.24,-1.10, 0.48,...    5,54,131,...       P,Y,X,...
#> 100       100    2.13,0.19,1.29,...   15,38,461,...       L,B,E,...
#>                              ir         r     invar
#>                   <IRangesList> <RleList> <numeric>
#> 1     30-39,165-174,295-304,... 1,1,5,...         2
#> 2     90-99,267-276,411-420,... 5,2,3,...         3
#> 3      4-13,112-121,449-458,... 5,4,4,...         4
#> 4   166-175,291-300,293-302,... 3,1,1,...         5
#> 5      35-44,93-102,106-115,... 5,1,2,...         6
#> ...                         ...       ...       ...
#> 96  233-242,382-391,385-394,... 4,5,1,...        97
#> 97         9-18,14-23,76-85,... 3,4,5,...        98
#> 98  116-125,135-144,193-202,... 4,1,4,...        99
#> 99       5-14,54-63,131-140,... 2,5,2,...       100
#> 100     15-24,38-47,461-470,... 5,2,2,...       101

## With a tally of the number of members in each group
reduceDataFrame(df, df$k, count = TRUE)
#> DataFrame with 100 rows and 8 columns
#>             k                     x               y               z
#>     <integer>         <NumericList>   <IntegerList> <CharacterList>
#> 1           1 -0.69,-0.15,-1.90,...  30,165,295,...       P,I,M,...
#> 2           2  0.33, 1.16,-0.88,...  90,267,411,...       P,S,W,...
#> 3           3 -0.77,-0.14, 0.65,...   4,112,449,...       V,M,T,...
#> 4           4    0.73,0.52,1.67,... 166,291,293,...       E,U,H,...
#> 5           5 -0.47,-0.38,-0.10,...   35,93,106,...       K,C,P,...
#> ...       ...                   ...             ...             ...
#> 96         96 -0.32,-0.63, 1.09,... 233,382,385,...       U,U,B,...
#> 97         97  0.79, 0.43,-0.66,...     9,14,76,...       I,T,Z,...
#> 98         98  0.83,-0.53, 0.25,... 116,135,193,...       A,C,S,...
#> 99         99 -1.24,-1.10, 0.48,...    5,54,131,...       P,Y,X,...
#> 100       100    2.13,0.19,1.29,...   15,38,461,...       L,B,E,...
#>                              ir         r     invar        .n
#>                   <IRangesList> <RleList> <numeric> <integer>
#> 1     30-39,165-174,295-304,... 1,1,5,...         2        13
#> 2     90-99,267-276,411-420,... 5,2,3,...         3         9
#> 3      4-13,112-121,449-458,... 5,4,4,...         4         9
#> 4   166-175,291-300,293-302,... 3,1,1,...         5        11
#> 5      35-44,93-102,106-115,... 5,1,2,...         6         8
#> ...                         ...       ...       ...       ...
#> 96  233-242,382-391,385-394,... 4,5,1,...        97         6
#> 97         9-18,14-23,76-85,... 3,4,5,...        98        12
#> 98  116-125,135-144,193-202,... 4,1,4,...        99        10
#> 99       5-14,54-63,131-140,... 2,5,2,...       100        16
#> 100     15-24,38-47,461-470,... 5,2,2,...       101        11

## Much faster, but more crowded result
df3 <- reduceDataFrame(df, df$k, simplify = FALSE)
df3
#> DataFrame with 100 rows and 7 columns
#>                   k                     x               y               z
#>       <IntegerList>         <NumericList>   <IntegerList> <CharacterList>
#> 1         1,1,1,... -0.69,-0.15,-1.90,...  30,165,295,...       P,I,M,...
#> 2         2,2,2,...  0.33, 1.16,-0.88,...  90,267,411,...       P,S,W,...
#> 3         3,3,3,... -0.77,-0.14, 0.65,...   4,112,449,...       V,M,T,...
#> 4         4,4,4,...    0.73,0.52,1.67,... 166,291,293,...       E,U,H,...
#> 5         5,5,5,... -0.47,-0.38,-0.10,...   35,93,106,...       K,C,P,...
#> ...             ...                   ...             ...             ...
#> 96     96,96,96,... -0.32,-0.63, 1.09,... 233,382,385,...       U,U,B,...
#> 97     97,97,97,...  0.79, 0.43,-0.66,...     9,14,76,...       I,T,Z,...
#> 98     98,98,98,...  0.83,-0.53, 0.25,... 116,135,193,...       A,C,S,...
#> 99     99,99,99,... -1.24,-1.10, 0.48,...    5,54,131,...       P,Y,X,...
#> 100 100,100,100,...    2.13,0.19,1.29,...   15,38,461,...       L,B,E,...
#>                              ir         r           invar
#>                   <IRangesList> <RleList>   <NumericList>
#> 1     30-39,165-174,295-304,... 1,1,5,...       2,2,2,...
#> 2     90-99,267-276,411-420,... 5,2,3,...       3,3,3,...
#> 3      4-13,112-121,449-458,... 5,4,4,...       4,4,4,...
#> 4   166-175,291-300,293-302,... 3,1,1,...       5,5,5,...
#> 5      35-44,93-102,106-115,... 5,1,2,...       6,6,6,...
#> ...                         ...       ...             ...
#> 96  233-242,382-391,385-394,... 4,5,1,...    97,97,97,...
#> 97         9-18,14-23,76-85,... 3,4,5,...    98,98,98,...
#> 98  116-125,135-144,193-202,... 4,1,4,...    99,99,99,...
#> 99       5-14,54-63,131-140,... 2,5,2,... 100,100,100,...
#> 100     15-24,38-47,461-470,... 5,2,2,... 101,101,101,...

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
