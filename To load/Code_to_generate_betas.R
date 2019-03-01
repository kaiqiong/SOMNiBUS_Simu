
# Remove the patient ID and sample identifier from the data set.

all(colnames(dat.use.total)==rownames(my.pheno))

rownames(my.pheno) <- colnames(dat.use.total) <- paste0("Sample", 1:43)

my.pheno <- my.pheno[,-3]

save(pos, beta.0, beta.1, beta.2,dat.use.total,  file="BANK1betas.RData")
