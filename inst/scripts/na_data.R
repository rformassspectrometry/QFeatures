m <- matrix(1:12, ncol = 3)
colnames(m) <- LETTERS[1:3]
rownames(m) <- letters[1:4]
m[3, 2] <- m[4, 1] <- m[1, 1] <- NA
se_na <- SummarizedExperiment(assay = m,
                              rowData = DataFrame(X = rep(1:2, 2)))
ft_na <- Features(list(na = se_na),
                  colData = DataFrame(row.names = LETTERS[1:3]))

save(ft_na, file = "../../data/ft_na.rda", compress = "xz", compression_level = 9)
