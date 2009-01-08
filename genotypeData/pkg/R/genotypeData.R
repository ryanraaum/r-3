CLASS = "genotypeData"

# --------------------------/
# Class Definition
# --------------------------/

setClass(CLASS,
	representation(
		genotypes = "matrix",
		sample_groups = "factor",
		sample_sizes = "numeric",
		ploidy  = "numeric",
		marker_groups = "factor",
		extras = "data.frame"
	),
	prototype = prototype(
		genotypes = matrix(nr=0,nc=0),
		sample_groups = factor(),
		sample_sizes = numeric(0),
		ploidy = numeric(0),
		marker_groups = factor(),
		extras = data.frame()
	)
)

# --------------------------/
# Validation
# --------------------------/

.validate <- function(object) { 
    msg <- NULL
    # samples are in the rows in the genotypes matrix,
    # so all information pertaining to samples must have
    # the same number of entries as there are rows
    #
    # -- groups -- every sample must be assigned to a group
    #              (even if every sample is assigned to the same group)
	if (length(samples(object)) > length(sampleGroups(object)))
   		msg <- c(msg, "there are more samples than group assignments")
    else if (length(samples(object)) < length(sampleGroups(object)))
    	msg <- c(msg, "there are more group assignments than samples")
    # -- sampleSizes -- every sample must have a sample size
    if (length(samples(object)) > length(sampleSizes(object)))
    	msg <- c(msg, "there are more samples than sample size data")
    else if (length(samples(object)) < length(sampleSizes(object)))
    	msg <- c(msg, "there are more sample size data than samples")

    if (is.null(msg)) TRUE
    else msg
}

setValidity(CLASS, .validate)

# --------------------------/
# Constructor
# --------------------------/

genotypeData <- function(genotypes,
						 ploidy=2,
						 onerowperind=FALSE,
						 samples=NULL,
						 sampleGroups=NULL,
						 sampleSizes=NULL,
						 markers=NULL,
						 markerGroups=NULL,
						 extras=data.frame()) {
	
	if (! is.matrix(genotypes))
		genotypes <- as.matrix(genotypes)

	if (! is.numeric(ploidy))
		stop("argument 'ploidy' of wrong type (must be numeric)")

	if (length(ploidy) != 1)
		stop("argument 'ploidy' of wrong length (must be a single number)")

	if ((ploidy %% 1) != 0)
		stop("argument 'ploidy' must be a whole number")

	if (! is.logical(onerowperind))
		stop("argument 'onerowperind' of wrong type (must be logical)")

	nrows <- dim(genotypes)[1]
	ncols <- dim(genotypes)[2]

	if (is.null(sampleGroups))
		sampleGroups <- factor(rep(1, nrows))

	if (! is.factor(sampleGroups))
		sampleGroups <- factor(sampleGroups)

	if (is.null(sampleSizes))
		sampleSizes <- rep(1, nrows)

	if (is.null(markerGroups))
		markerGroups <- factor(rep(1, ncols))

	if (! is.factor(markerGroups))
		markerGroups <- factor(markerGroups)

	if (onerowperind) {

		# with one row per individual, the number of columns in the 
		# genotypes matrix must be evenly divisible by the ploidy
		if (ncols %% ploidy != 0)
			stop("the number of columns in the genotypes data is not compatible with the given ploidy")

		if (is.null(samples))
			samples <- seq(1, nrows)

		if (is.null(markers))
			markers <- rep(seq(1, ncols/ploidy), each=ploidy)

		if (length(markers) == ncols/ploidy)
			markers <- rep(markers, each=ploidy)

		if (length(markerGroups) == ncols/ploidy)
			markerGroups <- rep(markerGroups, each=ploidy)

	} else {

		# with multiple rows per individual, the number of rows in the 
		# genotypes matrix must be evenly divisible by the ploidy
		if (nrows %% ploidy != 0)
			stop("the number of rows in the genotypes data is not compatible with the given ploidy")

		if (is.null(samples))
			samples <- rep(seq(1, nrows/ploidy), each=ploidy)

		if (length(samples) == nrows/ploidy)
			samples <- rep(samples, each=ploidy)

		if (length(sampleGroups) == nrows/ploidy)
			sampleGroups <- rep(sampleGroups, each=ploidy)

		if (length(sampleSizes) == nrows/ploidy)
			sampleSizes <- rep(sampleSizes, each=ploidy)

		if (is.null(markers))
			markers <- seq(1,ncols)
	}

	if (nrows != length(samples))
		stop("the number of samples is not consistent with the genotype data matrix")

	if (length(samples) != length(sampleGroups))
		stop("the number of sample groups is not consistent with the number of samples")

	if (length(samples) != length(sampleSizes))
		stop("the number of sample sizes is not consistent with the number of samples")

	if (ncols != length(markers))
		stop("the number of markers is not consistent with the genotype data matrix")

	if (onerowperind) {

		# all of this is to ensure that each marker name and each marker group is 
		# repeated ploidy times
		# i.e. for ploidy=2 the marker list is something like c("M1", "M1", "M2", "M2")
		selectors <- data.frame(row.names = seq(1, ncols/ploidy))
		for (i in 1:ploidy)
			selectors <- cbind(selectors, seq(i, ncols, by=ploidy))
		last <- 1
		for (current in 2:ploidy) {
			if (! all(markers[selectors[,last]] == markers[selectors[,current]]))
				stop("markers vector provided is not compatible with genotypes matrix")
			if (! all(markerGroups[selectors[,last]] == markerGroups[selectors[,current]]))
				stop("marker groups vector provided is not compatible with genotypes matrix")
			last <- current
		}

	} else {

		# all of this is to ensure that each sample name is repeated ploidy times
		# i.e. for ploidy=2 the sample list is something like c("S1", "S1", "S2", "S2")
		selectors <- data.frame(row.names = seq(1, nrows/ploidy))
		for (i in 1:ploidy)
			selectors <- cbind(selectors, seq(i, nrows, by=ploidy))
		last <- 1
		for (current in 2:ploidy) {
			if (! all(samples[selectors[,last]] == samples[selectors[,current]]))
				stop("sample vector provided is not compatible with genotypes matrix")
			if (! all(sampleGroups[selectors[,last]] == sampleGroups[selectors[,current]]))
				stop("sampleGroups vector provided is not compatible with genotypes matrix")
			if (! all(sampleSizes[selectors[,last]] == sampleSizes[selectors[,current]]))
				stop("sampleSizes vector provided is not compatible with genotypes matrix")
			last <- current
		}

	}

	rownames(genotypes) <- samples
	colnames(genotypes) <- markers

	new("genotypeData",
		genotypes=genotypes,
		sample_groups=sampleGroups,
		sample_sizes=sampleSizes,
		marker_groups=markerGroups,
		extras=extras)
}

# --------------------------/
# Utility Functions
# --------------------------/

# transform genotypes data matrix from 
# multiple rows per individual to one row per individual
# (using matrix algebra magic)
.multi2one <- function(genotypes, ploidy) {
	nrows <- dim(genotypes)[1]
	nsamples <- nrows/ploidy
	A <- t(kronecker(diag(1,nsamples),as.vector(diag(1,ploidy)))) 
	B <- kronecker(genotypes,diag(1,ploidy))
	A %*% B
}

# and vice versa
.one2multi <- function(genotypes, ploidy) {
	nsamples <- dim(genotypes)[1]
	A <- kronecker(diag(1,nsamples+1),as.vector(diag(1,ploidy)))
	B <- kronecker(genotypes,diag(1,ploidy))
	B %*% A
}


# --------------------------/
# Accessors - Getters
# --------------------------/

bindMethod("genotypes",   	CLASS, function(object) object@genotypes)
bindMethod("samples",     	CLASS, function(object) rownames(object@genotypes))
bindMethod("sampleGroups",	CLASS, function(object) object@sample_groups)
bindMethod("sampleSizes", 	CLASS, function(object) object@sample_sizes)
bindMethod("ploidy",      	CLASS, function(object) object@ploidy)
bindMethod("markers",     	CLASS, function(object) colnames(object@genotypes))
bindMethod("markerGroups",	CLASS, function(object) object@marker_groups)
bindMethod("extras",		CLASS, function(object) object@extras)

# --------------------------/
# Accessors - Setters
# --------------------------/

bindReplaceMethod("genotypes<-", CLASS, function(object, value) {
  object@genotypes <- value
  object
})

bindReplaceMethod("samples<-", CLASS, function(object, value) {
  rownames(object@genotypes) <- value
  object
})

bindReplaceMethod("sampleGroups<-", CLASS, function(object, value) {
  object@sample_groups <- value
  object
})

bindReplaceMethod("sampleSizes<-", CLASS, function(object, value) {
  object@sample_sizes <- value
  object
})

bindReplaceMethod("ploidy<-", CLASS, function(object, value) {
  object@ploidy <- value
  object
})

bindReplaceMethod("markers<-", CLASS, function(object, value) {
  colnames(object@genotypes) <- value
  object
})

bindReplaceMethod("markerGroups<-", CLASS, function(object, value) {
  object@marker_groups <- value
  object
})

bindReplaceMethod("extras<-", CLASS, function(object, value) {
  object@extras <- value
  object
})

# --------------------------/
# Extract Subsets
# --------------------------/

bindMethod("extract", CLASS, function(object, grps) extr(object, grps))

"markerColumns" <-
function (gd, mrkrs)
{
	positions <- unlist(mapply(rep, markers(gd), ploidy(gd)), use.names=FALSE)
	apply_func <- function (x) which(x == positions)
	return(as.vector(sapply(mrkrs, apply_func, simplify=TRUE, USE.NAMES=FALSE)))
}

"extr" <-
function (gd, grps=NULL, mrkrs=NULL)
{
	select_g <- NULL
	if (is.null(grps)) {
		select_g <- rep(TRUE, length(samples(gd)))
	}
	else {
		select_g <- rep(FALSE, length(samples(gd)))
		for (i in 1:length(grps)) {
			select_g <- select_g | (groups(gd) == grps[i])
		}
	}

	select_m  <- NULL
	select_mc <- NULL
	if (is.null(mrkrs)) {
		select_m  <- 1:length(markers(gd))
		select_mc <- markerColumns(gd, markers(gd))
	}
	else {
		gd_markers <- markers(gd)
		apply_func <- function (x) match(x, gd_markers)
		select_m  <- sapply(markers(gd), apply_func, simplify=TRUE, USE.NAMES=FALSE)
		select_mc <- markerColumns(gd, markers(gd)[select_m])
	}
	
	ngd <- new( CLASS,
				ploidy=ploidy(gd)[select_m],
				genotypes=genotypes(gd)[select_g,select_mc],
				sample_sizes=sampleSizes(gd)[select_g],
				groups=factor(groups(gd)[select_g]),
				markers=markers(gd)[select_m] )
				
	return(ngd)
}

