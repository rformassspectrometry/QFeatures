m <- matrix(1:12, ncol = 3)
colnames(m) <- LETTERS[1:3]
rownames(m) <- letters[1:4]
m[3, 2] <- m[4, 1] <- m[1, 1] <- NA
se_na <- SummarizedExperiment(assay = m,
                              rowData = DataFrame(X = rep(1:2, 2),
                                                  Y = rep(LETTERS[1:2], 2)))
ft_na <- QFeatures(list(na = se_na),
                  colData = DataFrame(row.names = LETTERS[1:3]))

save(ft_na, file = "../../data/ft_na.rda", compress = "xz", compression_level = 9)


library("MSnbase")
library("pRolocdata")
data(naset, package = "pRolocdata")

se_na2 <- as(naset, "SummarizedExperiment")
ft_na2 <- QFeatures(list(x = se_na2), colData = colData(se_na2))

save(se_na2, file = "../../data/se_na2.rda")
## save(ft_na2, file = "../../data/ft_na2.rda")
