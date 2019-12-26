## Same data as in the example
k <- sample(100, 1e3, replace = TRUE)
df <- DataFrame(k = k,
                x = round(rnorm(length(k)), 2),
                y = seq_len(length(k)),
                z = sample(LETTERS, length(k), replace = TRUE),
                ir = IRanges(seq_along(k), width = 10),
                r = Rle(sample(5, length(k), replace = TRUE)),
                invar = k + 1)


test_that("reduce DataFrame", {
    df2 <- reduceDataFrame(df, df$k)
    expect_identical(dim(df2), c(100L, ncol(df)))
    expect_identical(unname(df2[["k"]]), 1:100)

    expect_identical(unlist(df2[1, "y"], use.names = FALSE),
                     df[df$k == 1, "y"])
    expect_identical(unlist(df2[2, "x"], use.names = FALSE),
                     df[df$k == 2, "x"])
    expect_identical(unlist(df2[100, "z"], use.names = FALSE),
                     df[df$k == 100, "z"])
    expect_identical(unname(df2[100, "invar"]),
                     unique(df[df$k == 100, "invar"]))    
    df3 <- reduceDataFrame(df, df$k, count = TRUE)
    expect_identical(dim(df3), c(100L, ncol(df) + 1L))
    df4 <- reduceDataFrame(df, df$k, simplify = FALSE)
    expect_identical(dim(df4), c(100L, ncol(df)))
    expect_identical(df2[, -c(1, 7)], df4[, -c(1, 7)])
    df5 <- reduceDataFrame(df, df$k, drop = TRUE)
    expect_identical(dim(df5), c(100L, 2L))
})


test_that("expand DataFrame", {
    df2 <- reduceDataFrame(df, df$k, simplify = FALSE, count = FALSE)
    df0 <- expandDataFrame(df2)
    expect_identical(dim(df0), dim(df))
    expect_identical(names(df0), names(df))
    df0_y <- df0[order(df0$y), ]
    df_y <- df[order(df$y), ]
    rownames(df0_y) <- NULL
    expect_identical(df0_y, df_y)
})


test_that("expand DataFrame (with k)", {
    df2 <- reduceDataFrame(df, df$k, simplify = FALSE, count = FALSE)
    df0 <- expandDataFrame(df2, df$k)
    expect_identical(dim(df0), dim(df))
    expect_identical(names(df0), names(df))
    df0_y <- df0[order(df0$y), ]
    df_y <- df[order(df$y), ]
    rownames(df0_y) <- NULL
    expect_identical(df0_y, df_y)
})
