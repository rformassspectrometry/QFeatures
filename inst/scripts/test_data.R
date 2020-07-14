library("QFeatures")

set.seet(1)

## ---------------------------------------------------------
## Creating psms as <double> rather than <integer>
## ---------------------------------------------------------
psms <- matrix(seq(1, 20, by = 1.0), ncol = 2)
colnames(psms) <- paste0("S", 1:2)
rowdata <- DataFrame(Sequence = c("SYGFNAAR", "SYGFNAAR", "SYGFNAAR", "ELGNDAYK",
                                  "ELGNDAYK", "ELGNDAYK", "IAEESNFPFIK",
                                  "IAEESNFPFIK", "IAEESNFPFIK", "IAEESNFPFIK"),
                     Protein = c("ProtA", "ProtA", "ProtA", "ProtA", "ProtA",
                                 "ProtA", "ProtB", "ProtB", "ProtB", "ProtB"),
                     Var = 1:10,
                     location = c(rep("Mitochondrion", 6), rep("unknown", 4)),
                     pval = round(runif(10, max = 0.1, min = 0), 3))
rownames(rowdata) <- rownames(psms) <- paste0("PSM", 1:10)
coldata <- DataFrame(Group = 1:2)
rownames(coldata) <- colnames(psms)


psms <- SummarizedExperiment(psms, rowData = rowdata)

feat1 <- QFeatures(list(psms = psms), colData = coldata)

save(feat1, file = "../../data/feat1.rda", compress = "xz", compression_level = 9)


## ---------------------------------------------------------
## Test data for joining assays
## ---------------------------------------------------------

m1 <- matrix(1:40, ncol = 4) + 0.0
m2 <- matrix(1:16, ncol = 4) + 0.0
m3 <- matrix(1:28, ncol = 4) + 0.0
colnames(m1) <- paste0("S", 1:4)
colnames(m2) <- paste0("S", 5:8)
colnames(m3) <- paste0("S", 9:12)
rownames(m1) <- letters[1:10]
rownames(m2) <- letters[8:11]
rownames(m3) <- letters[c(1, 2, 10:14)]

df1 <- DataFrame(Prot = paste0("P", rownames(m1)),
                 x = rnorm(1:10),
                 row.names = rownames(m1))
cd1 <- DataFrame(Var1 = rnorm(4),
                 Var2 = LETTERS[1:4],
                 row.names = colnames(m1))


df2 <- DataFrame(Prot = paste0("P", rownames(m2)),
                 x = rnorm(1:4),
                 y = rnorm(1:4),
                 row.names = rownames(m2))
cd2 <- DataFrame(Var1 = rnorm(4),
                 Var2 = LETTERS[5:8],
                 row.names = colnames(m2))

df3 <- DataFrame(Prot = paste0("P", rownames(m3)),
                 x = rnorm(1:7),
                 y = rnorm(1:7),
                 row.names = rownames(m3))
cd3 <- DataFrame(Var1 = rnorm(4),
                 row.names = colnames(m3))


se1 <- SummarizedExperiment(m1, df1, colData = cd1)
se2 <- SummarizedExperiment(m2, df2, colData = cd2)
se3 <- SummarizedExperiment(m3, df3, colData = cd3)

el <- list(assay1 = se1, assay2 = se2, assay3 = se3)
feat2 <- QFeatures(el)

save(feat2, file = "../../data/feat2.rda", compress = "xz", compression_level = 9)
