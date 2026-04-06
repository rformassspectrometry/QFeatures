data(feat3)

se1 <- feat3[[1]]
rownames(se1) <- 1:nrow(se1)
rowData(se1)$charge <- c(1, 1, 2, 1, 2, 3, 1)

se2 <- feat3[[2]]
rownames(se2) <- 1:nrow(se2)
rowData(se2)$charge <- c(1, 1, 2, 3, 1, 1, 2, 3)

rowData(se1)$pval <- rowData(se2)$pval <- NULL
rowData(se1)$Var <- rowData(se2)$Var <- NULL
rowData(se1)$location <- rowData(se2)$location <- NULL

feat4 <- QFeatures(list(PSM1 = se1, PSM2 = se2))

stopifnot(validObject(feat4))

save(feat4, file = "../../data/feat4.rda")