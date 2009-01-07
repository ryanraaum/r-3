"readStructureData" <-
function (filename, 
		  header=TRUE, 
		  popdata=TRUE, 
		  popflag=TRUE, 
		  phenotype=FALSE, 
		  nextra=0, 
		  onerowperind=FALSE,
		  na.s="-9")
{	
	## What columns are present? Set up the first few column labels using flags
	pre <- c("Label")
	xtras <- c()
	if (popdata)   { pre <- c(pre, "PopData") }
	if (popflag)   { pre <- c(pre, "PopFlag"); xtras <- c(xtras, "PopFlag") }
	if (phenotype) { pre <- c(pre, "Phenotype"); xtras <- c(xtras, "Phenotypes") }
	if (nextra >= 1) {
		for (i in 1:nextra) { 
			name = paste("Extra", i, sep="")
			pre <- c(pre, name) 
			xtras <- c(xtras, name) 
		}
	}
	
	## When the main data are read in, make sure to skip the header if there is one
	nskip <- 0
	if (header) { nskip <- 1 }
	dat <- read.table(filename, skip=nskip, na.strings=na.s)

	genotypes <- as.matrix(dat[,(length(pre)+1):ncol(dat)])

	## If there is a header with the marker names, read that in 
	## otherwise just create some generic "V" marker column names
	if (header) {
		cnames <- as.vector(read.table(filename, nrows=1, colClasses="character"), 
							mode="character")
	} else {
		# if no marker names are given, we assume haploidy
		cnames <- paste("V", 1:(dim(genotypes)[2]), sep="")
	}

	rownames(genotypes) <- dat[,1]
	colnames(genotypes) <- cnames
	names(dat) <- c(pre, cnames)

	# use popdata as sample_groups, or lacking that, put all samples in one group
	if (popdata) {
		sample_groups <- factor(dat[,"PopData"])
	} else {
		sample_groups <- factor(rep(1, dim(genotypes)[1]))
	}
	
	# to determine ploidy, either sample names or marker names are repeated
	if (onerowperind) {
		# if one row per individual, then marker names should be repeated
		mn <- colnames(genotypes) # mn = marker names
		# we assume all markers have the same ploidy,
		# so check any marker - here the first.
		ploidy <- length(mn[mn==mn[1]])	
	} else {
		# if ploidy > 1 is indicated by multiple rows per individual, then
		sn <- colnames(genotypes) # sn = sample names
		# same assumptions as just above
		ploidy <- length(sn[sn==sn[1]])
	}

	gd = new("genotypeData", 
			 genotypes=genotypes, 
			 extras=dat[,xtras],
			 ploidy=ploidy,
			 sample_groups=sample_groups,
			 sample_sizes=rep(1,dim(genotypes)[1]),
			 marker_groups=factor(rep(1,dim(genotypes)[2])))

	return(gd)
}

