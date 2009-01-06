### --- Test setup ---

a <- NULL

.setUp <- function () a <<- new("genotypeData")

### --- Test functions ---
 
test.default_values <- function()
{
  checkEquals(matrix(nr=0,nc=0), genotypes(a))
  checkEquals(NULL, samples(a))
  checkEquals(factor(), sampleGroups(a))
  checkEquals(numeric(0), sampleSizes(a))
  checkEquals(numeric(0), ploidy(a))
  checkEquals(character(0), markers(a))
  checkEquals(factor(), markerGroups(a))
}

test.assignment <- function()
{
  genotypes(a) <- matrix(c(1,1,2,2), nrow=1)
  checkEquals(1, nrow(genotypes(a)))
  checkEquals(4, ncol(genotypes(a)))

  samples(a) <- "S1"
  checkEquals("S1", samples(a))

  sampleGroups(a) <- factor(c(1,1))
  checkEquals(factor(c(1,1)), sampleGroups(a))

  sampleSizes(a) <- 1
  checkEquals(1, sampleSizes(a))

  ploidy(a) <- c(2,2)
  checkEquals(c(2,2), ploidy(a))

  markers(a) <- c("M1", "M2")
  checkEquals(c("M1", "M2"), markers(a))

  markerGroups(a) <- factor(c(1,1))
  checkEquals(factor(c(1,1)), markerGroups(a))
}

test.validation <- function()
{
    # the heart of the genotypeData class is a matrix of genotype data
    genotypes <- matrix(c(1,1,2,2), nr=1)
    # the sample names are the rownames of the genotypes matrix
    rownames(genotypes) <- "S1"

    # this should throw an exception because 
    # 1. groups and sampleSizes
    #    need to have information to match the number of samples
    # 2. ploidy and markers
    #    need to match up to the genotype data columns
    checkException(new("genotypeData", genotypes=genotypes))

    # first create some group and sample size data to match the genotype data
    sample_groups <- factor(c(1))
    sample_sizes <- 1

    # next create some ploidy and marker data to match the genotype data
    ploidy <- c(2,2)
    markers <- c("M1", "M2")
    marker_groups <- factor(c(1,1))

    # now this should pass
    new("genotypeData", genotypes=genotypes, 
                        sample_groups=sample_groups,
                        sample_sizes=sample_sizes,
                        ploidy=ploidy,
                        markers=markers,
						marker_groups=marker_groups)
}

