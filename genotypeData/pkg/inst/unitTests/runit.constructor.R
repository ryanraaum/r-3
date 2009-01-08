### --- Test setup ---

valid_diploid_multi_row_genotypes <- matrix(c(1,1,2,2,1,1,1,2), nr=4, byrow=T)

### --- Test functions ---
 
test.bad_input <- function()
{
	# requires valid genotype data as the first argument
	checkException(genotypeData())

	# genotype data cannot be NULL
	checkException(genotypeData(NULL))

	# abbreviate valid genotype matrix for convenience
	g <- valid_diploid_multi_row_genotypes

	# ploidy must be numeric
	checkException(genotypeData(g, ploidy="cat"))

	# ploidy must be a single number
	checkException(genotypeData(g, ploidy=c(2,2)))

	# ploidy must be a whole number
	checkException(genotypeData(g, ploidy=2.4))

	# onerowperind must be logical
	checkException(genotypeData(g, onerowperind=4))
	checkException(genotypeData(g, onerowperind="cat"))
	checkException(genotypeData(g, onerowperind=NULL))

	# for haploid input, sample names (if provided) must match the 
	# number of rows in the genotype data input
	checkException(genotypeData(g, ploidy=1, samples=1:5))
	checkException(genotypeData(g, ploidy=1, samples=1:3))

	# for diploid data, samples vector must have either:
	# 1. same number of entries as rows in genotype matrix, but doubled names
	# 2. half the number of entries as rows in genotype matrix
	checkException(genotypeData(g, samples=c(1,2,2,2)))
	checkException(genotypeData(g, samples=c(1,2,3,4)))

	# for haploid input, marker names (if provided) must match the 
	# number of columns in the genotype data input
	checkException(genotypeData(g, ploidy=1, markers=1:3))
	checkException(genotypeData(g, ploidy=1, markers=1))
}

test.good_input <- function()
{
	# abbreviate valid genotype matrix for convenience
	g <- valid_diploid_multi_row_genotypes

	# valid genotype matrix
	checkEquals("genotypeData", class(genotypeData(g))[[1]])

	# for haploid data, samples vector must have same
	# number of entries as rows in genotype matrix
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, samples=1:4))[[1]])

	# for diploid data, samples vector must have either:
	# 1. same number of entries as rows in genotype matrix, but doubled names
	# 2. half the number of entries as rows in genotype matrix
	checkEquals("genotypeData", class(genotypeData(g, samples=c(1,1,2,2)))[[1]])
	checkEquals("genotypeData", class(genotypeData(g, samples=c(1,2)))[[1]])

}
