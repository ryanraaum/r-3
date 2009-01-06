### --- Test setup ---

library(RArlequin)

### --- Test functions ---
 
test.functions <- function()
{
  checkTrue(is.function(readArlequinData))
  checkEquals(8, conv(2,4))
  checkEquals(2, testit("monkey"))
  checkEquals(2, testit2(2))
}

