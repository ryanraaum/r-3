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

	# should fail if the format of the matrix is not as expected
	checkException(genotypeData(g, ploidy=2, samples=c(1,2), onerowperind=T))
}

test.good_input.haploid <- function()
{
	# haploid is ploidy=1

	# abbreviate valid genotype matrix for convenience
	g <- valid_diploid_multi_row_genotypes

	# haploid with everything else default
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1))[[1]])
	# for haploid data, samples vector must have same
	# number of entries as rows in genotype matrix
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, samples=1:4))[[1]])
	# should still work with valid markers vector
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, samples=1:4, markers=1:2))[[1]])
	# should still work with only markers vector
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, markers=1:2))[[1]])
	# now with sampleGroups
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, sampleGroups=factor(c(1,1,2,2))))[[1]])
	# should silently deal with non-factor sample groups
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, sampleGroups=c(1,1,2,2)))[[1]])
	# now with markerGroups
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, markerGroups=factor(c(1,1))))[[1]])
	# should silently deal with non-factor sample groups
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, markerGroups=c(1,1)))[[1]])
	# now with sampleSizes
	checkEquals("genotypeData", class(genotypeData(g, ploidy=1, sampleSizes=c(1,1,1,1)))[[1]])
	# now with everything
	checkEquals("genotypeData", class(genotypeData(g, 
												   ploidy=1, 
												   samples=1:4,
												   sampleGroups=factor(c(1,1,2,2)),
												   sampleSizes=rep(1,4),
												   markers=1:2,
												   markerGroups=factor(c(1,1))))[[1]])
	# finally, for haploid, shouldn't actually matter if one sets onerowperind or not
	checkEquals("genotypeData", class(genotypeData(g, 
												   onerowperind=T,
												   ploidy=1, 
												   samples=1:4,
												   sampleGroups=factor(c(1,1,2,2)),
												   sampleSizes=rep(1,4),
												   markers=1:2,
												   markerGroups=factor(c(1,1))))[[1]])

}

test.good_input.diploid <- function()
{
	# diploid is ploidy=2
	# which is the default ploidy value for the constructor

	# abbreviate valid genotype matrix for convenience
	g <- valid_diploid_multi_row_genotypes
	# for these tests, this genotype matrix has 2 diploid samples for 2 markers

	# FIRST - with multiple rows per individual (default)

	# valid genotype matrix
	checkEquals("genotypeData", class(genotypeData(g))[[1]])
	# for diploid data, samples vector must have either 
	# 1. include each sample name once 
	checkEquals("genotypeData", class(genotypeData(g, samples=1:2))[[1]])
	# 2. include each sample name 2 times (for each row)
	checkEquals("genotypeData", class(genotypeData(g, samples=c(1,1,2,2)))[[1]])
	# should still work with valid markers vector
	checkEquals("genotypeData", class(genotypeData(g, samples=1:2, markers=1:2))[[1]])
	# should still work with only markers vector
	checkEquals("genotypeData", class(genotypeData(g, markers=1:2))[[1]])
	# now with sampleGroups
	checkEquals("genotypeData", class(genotypeData(g, sampleGroups=factor(c(1,1,2,2))))[[1]])
	# should silently deal with non-factor sample groups
	checkEquals("genotypeData", class(genotypeData(g, sampleGroups=c(1,1,2,2)))[[1]])
	# now with markerGroups
	checkEquals("genotypeData", class(genotypeData(g, markerGroups=factor(c(1,1))))[[1]])
	# should silently deal with non-factor sample groups
	checkEquals("genotypeData", class(genotypeData(g, markerGroups=c(1,1)))[[1]])
	# now with sampleSizes
	checkEquals("genotypeData", class(genotypeData(g, sampleSizes=c(1,1,1,1)))[[1]])
	# now with everything
	checkEquals("genotypeData", class(genotypeData(g, 
												   samples=c(1,1,2,2),
												   sampleGroups=factor(c(1,1,2,2)),
												   sampleSizes=rep(1,4),
												   markers=1:2,
												   markerGroups=factor(c(1,1))))[[1]])


	# NEXT - with a single row per individual (onerowperind=TRUE)
}
