"readStructureData" <-
function (filename, header=TRUE, popdata=TRUE, popflag=TRUE, phenotype=FALSE, nextra=0, na.s="-9")
{	
	## What columns are present? Set up the first few column labels using flags
	pre <- c("Label")
	if (popdata)   { pre <- c(pre, "PopData") }
	if (popflag)   { pre <- c(pre, "PopFlag") }
	if (phenotype) { pre <- c(pre, "Phenotype") }
	pre <- c(pre, "NumMarkers")
	if (nextra >= 1) {
		for (i in 1:nextra) { 
			pre <- c(pre, paste("Extra", i, sep="")) 
		}
	}
	
	## When the main data are read in, make sure to skip the header if there is one
	nskip <- 0
	if (header) { nskip <- 1 }
	dat <- read.table(filename, skip=nskip, na.strings=na.s)
	
	## extra -1 because I've added the NumMarkers name
	lpre = length(pre) - 1
	nmarkers = length(dat[1,]) - lpre
	dat <- cbind(dat[1:lpre], nmarkers, dat[(lpre+1):length(dat[1,])])

	## If there is a header with the marker names, read that in and use it,
	## otherwise just create some generic "V" marker column names
	if (header) {
		cnames <- as.vector(read.table(filename, nrows=1, colClasses="character"), 
							mode="character")
		names(dat) <- c(pre, cnames)
	} else {
		names(dat) <- c(pre, paste("V", 1:nmarkers, sep=""))
	}
	
	return(dat)
}

