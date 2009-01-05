### --- Test setup ---
 
a <- new("genotypeData")
 
### --- Test functions ---
 
test.instantiation <- function()
{
  checkEquals(character(0), markers(a))
}


