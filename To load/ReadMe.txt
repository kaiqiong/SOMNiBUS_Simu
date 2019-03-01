This "To load" folder contains the functions and RData that will be required to "load" or "source" in order to excute the R scripts under other folders

------------------
Functions-V12.R
------------------
The functions inside file "Functions-V12.R" have been well documented in R package SOMNiBUS. See function documentation inside the package for a detailed overview.


------------------
BANK1betas.RData
-------------------
The objects stored in these BANK1betas.RData includes

1. dat.use.total: (A matrix with 123 rows and 43 columns)stores the read-depth in our real data example. Each row is a CpG site and Each columns is sample. this information will be used to generate the read-depth in our simulated data set.

2. my.pheno: (a data frame with 43 rows and 2 columns). The RA stauts  and cell type status for each sample. The proportions of RA and Tcell samples in our illustrative dataset were used to simulate covarites in our simulation experiments.

3. beta.0: the estimated beta.0(t) of each CpG sites at position t. (Intercept). Reseults obtained by applying SOMNiBUS to our real-life data example.  (detailed steps in SOMNiBUS vignette )

4. beta.1: the estimated beta.1(t) of each CpG sites at position t. (Effect of RA). Reseults obtained by applying SOMNiBUS to our real-life data example.  (detailed steps in SOMNiBUS vignette )

5. beta.2: the estimated beta.2(t) of each CpG sites at position t. (Effect of cell type). Reseults obtained by applying SOMNiBUS to our real-life data example.  (detailed steps in SOMNiBUS vignette )

6. pos: the genomic positions for the 123 CpG sites under investigation. The orders coinside the order in object betas.

e.g. plot(pos, beta.1) gives the smooth curve for the effect of RA

The shapes stored in beta.0, beta.1, and beta.2 were used to generate our simulated data
 