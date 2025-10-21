# Unfold a data frame

A data frame is said to be *folded* when some cells contain multiple
elements. These are often encode as a semi-colon separated character ,
such as `"a;b"`. This function will transform the data frame to that
`"a"` and `"b"` are split and recorded across two lines.

The simple example below illustrates a trivial case, where the table
below

|     |     |
|-----|-----|
| X   | Y   |
| 1   | a;b |
| 2   | c   |

is unfolded based on the Y variable and becomes

|     |     |
|-----|-----|
| X   | Y   |
| 1   | a   |
| 1   | b   |
| 2   | c   |

where the value 1 of variable X is now duplicated.

If there is a second variable that follows the same pattern as the one
used to unfold the table, it also gets unfolded.

|     |     |     |
|-----|-----|-----|
| X   | Y   | Z   |
| 1   | a;b | x;y |
| 2   | c   | z   |

becomes

|     |     |     |
|-----|-----|-----|
| X   | Y   | Z   |
| 1   | a   | x   |
| 1   | b   | y   |
| 2   | c   | z   |

because it is implied that the element in "a;b" are match to "x;y" by
their respective indices. Note in the above example, unfolding by Y or Z
produces the same result.

However, the following table unfolded by Y

|     |     |     |
|-----|-----|-----|
| X   | Y   | Z   |
| 1   | a;b | x;y |
| 2   | c   | x;y |

produces

|     |     |     |
|-----|-----|-----|
| X   | Y   | Z   |
| 1   | a   | x;y |
| 1   | b   | x;y |
| 2   | c   | x;y |

because "c" and "x;y" along the second row don't match. In this case,
unfolding by Z would produce a different result. These examples are also
illustrated below.

Note that there is no `foldDataFrame()` function. See
[`reduceDataFrame()`](https://rformassspectrometry.github.io/QFeatures/reference/reduceDataFrame.md)
and
[`expandDataFrame()`](https://rformassspectrometry.github.io/QFeatures/reference/reduceDataFrame.md)
to flexibly encode and handle vectors of length \> 1 within cells.

## Usage

``` r
unfoldDataFrame(x, k, split = ";")
```

## Arguments

- x:

  A `DataFrame` or `data.frame` to be unfolded.

- k:

  `character(1)` referring to a character variable in `x`, that will be
  used to unfold `x`.

- split:

  `character(1)` passed to
  [`strsplit()`](https://rdrr.io/r/base/strsplit.html) to split
  `x[[k]]`.

## Value

A new object unfolded object of class `class(x)` with numbers of rows
\>= `nrow(x)` and columns identical to `x`.

## Author

Laurent Gatto

## Examples

``` r

(x0 <- DataFrame(X = 1:2, Y = c("a;b", "c")))
#> DataFrame with 2 rows and 2 columns
#>           X           Y
#>   <integer> <character>
#> 1         1         a;b
#> 2         2           c
unfoldDataFrame(x0, "Y")
#> DataFrame with 3 rows and 2 columns
#>           X           Y
#>   <integer> <character>
#> 1         1           a
#> 2         1           b
#> 3         2           c

(x1 <- DataFrame(X = 1:2, Y = c("a;b", "c"), Z = c("x;y", "z")))
#> DataFrame with 2 rows and 3 columns
#>           X           Y           Z
#>   <integer> <character> <character>
#> 1         1         a;b         x;y
#> 2         2           c           z
unfoldDataFrame(x1, "Y")
#> DataFrame with 3 rows and 3 columns
#>           X           Y           Z
#>   <integer> <character> <character>
#> 1         1           a           x
#> 2         1           b           y
#> 3         2           c           z
unfoldDataFrame(x1, "Z") ## same
#> DataFrame with 3 rows and 3 columns
#>           X           Y           Z
#>   <integer> <character> <character>
#> 1         1           a           x
#> 2         1           b           y
#> 3         2           c           z

(x2 <- DataFrame(X = 1:2, Y = c("a;b", "c"), Z = c("x;y", "x;y")))
#> DataFrame with 2 rows and 3 columns
#>           X           Y           Z
#>   <integer> <character> <character>
#> 1         1         a;b         x;y
#> 2         2           c         x;y
unfoldDataFrame(x2, "Y")
#> DataFrame with 3 rows and 3 columns
#>           X           Y           Z
#>   <integer> <character> <character>
#> 1         1           a         x;y
#> 2         1           b         x;y
#> 3         2           c         x;y
unfoldDataFrame(x2, "Z") ## different
#> DataFrame with 4 rows and 3 columns
#>           X           Y           Z
#>   <integer> <character> <character>
#> 1         1         a;b           x
#> 2         1         a;b           y
#> 3         2           c           x
#> 4         2           c           y
```
