CLASS = "genotypeData"

# --------------------------/
# Class Definition
# --------------------------/

setClass(CLASS,
	representation(
		genotypes = "matrix",
		groups = "factor",
		sample_sizes = "numeric",
		ploidy  = "numeric",
		markers = "character",
        phased = "logical",
		description = "character",
		notes = "character"		
	),
	prototype = prototype(
		genotypes = matrix(nr=0,nc=0),
		groups = factor(),
		sample_sizes = numeric(0),
		ploidy = numeric(0),
		markers = character(0),
        phased = FALSE,
		description = "",
		notes = ""
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
    if (length(samples(object)) > length(groups(object)))
        msg <- c(msg, "there are more samples than group assignments")
    else if (length(samples(object)) < length(groups(object)))
        msg <- c(msg, "there are more group assignments than samples")
    # -- sampleSizes -- every sample must have a sample size
    if (length(samples(object)) > length(sampleSizes(object)))
        msg <- c(msg, "there are more samples than sample size data")
    else if (length(samples(object)) < length(groups(object)))
        msg <- c(msg, "there are more sample size data than samples")

    # marker data are in the columns of the genotypes matrix,
    # so the following must apply
    #
    # -- ploidy -- the sum of the ploidy vector must equal the number of 
    #              columns in the the genotypes matrix
    #              (each entry in the ploidy vector is the number of columns
    #               for the corresponding marker)
    if (sum(ploidy(object)) != ncol(genotypes(object)))
        msg <- c(msg, "ploidy values do not match up to the data")
    # -- markers -- there should be the same number of marker names as
    #               ploidy entries
    if (length(ploidy(object)) != length(markers(object)))
        msg <- c(msg, "markers and ploidy information do not match up")

    if (is.null(msg)) TRUE
    else msg
}

setValidity(CLASS, .validate)

# --------------------------/
# Accessors - Getters
# --------------------------/

bindMethod("ploidy",      CLASS, function(object) object@ploidy)

bindMethod("markers",     CLASS, function(object) object@markers)

bindMethod("genotypes",   CLASS, function(object) object@genotypes)

bindMethod("samples",     CLASS, function(object) rownames(object@genotypes))

bindMethod("groups",      CLASS, function(object) object@groups)

bindMethod("sampleSizes", CLASS, function(object) object@sample_sizes)

bindMethod("description", CLASS, function(object) object@description)

bindMethod("notes",       CLASS, function(object) object@notes)

# --------------------------/
# Accessors - Setters
# --------------------------/

bindReplaceMethod("ploidy<-", CLASS, function(object, value) {
  object@ploidy <- value
  object
})

bindReplaceMethod("markers<-", CLASS, function(object, value) {
  object@markers <- value
  object
})

bindReplaceMethod("genotypes<-", CLASS, function(object, value) {
  object@genotypes <- value
  object
})

bindReplaceMethod("samples<-", CLASS, function(object, value) {
  rownames(object@genotypes) <- value
  object
})

bindReplaceMethod("groups<-", CLASS, function(object, value) {
  object@groups <- value
  object
})

bindReplaceMethod("sampleSizes<-", CLASS, function(object, value) {
  object@sample_sizes <- value
  object
})

bindReplaceMethod("description<-", CLASS, function(object, value) {
  object@description <- value
  object
})

bindReplaceMethod("notes<-", CLASS, function(object, value) {
  object@notes <- value
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

