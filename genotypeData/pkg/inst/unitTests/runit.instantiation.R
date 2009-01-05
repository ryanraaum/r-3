### --- Test setup ---

a <- NULL

.setUp <- function () a <<- new("genotypeData")

### --- Test functions ---
 
test.default_values <- function()
{
  checkEquals(numeric(0), ploidy(a))
  checkEquals(character(0), markers(a))
  checkEquals(matrix(nr=0,nc=0), genotypes(a))
  checkEquals(factor(0), groups(a))
  checkEquals(numeric(0), sampleSizes(a))
  checkEquals("", description(a))
  checkEquals("", notes(a))
  checkEquals(NULL, samples(a))
}

test.assignment <- function()
{
  genotypes(a) <- matrix(c(1,1,2,2), nrow=1)
  checkEquals(1, nrow(genotypes(a)))
  checkEquals(4, ncol(genotypes(a)))

  ploidy(a) <- c(2,2)
  checkEquals(c(2,2), ploidy(a))

  markers(a) <- c("M1", "M2")
  checkEquals(c("M1", "M2"), markers(a))

  groups(a) <- factor(c(1,1))
  checkEquals(factor(c(1,1)), groups(a))

  sampleSizes(a) <- 1
  checkEquals(1, sampleSizes(a))

  samples(a) <- "S1"
  checkEquals("S1", samples(a))

  description(a) <- "blah"
  checkEquals("blah", description(a))

  notes(a) <- "blah"
  checkEquals("blah", notes(a))
}


