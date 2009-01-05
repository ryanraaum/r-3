CLASS = "genotypeData"

# --------------------------/
# Class Definition
# --------------------------/

setClass(CLASS,
	representation(
		ploidy  = "numeric",
		markers = "character",
		genotypes = "matrix",
		groups = "factor",
		sample_sizes = "numeric",
		description = "character",
		notes = "character"		
	),
	prototype = prototype(
		ploidy = numeric(0),
		markers = character(0),
		genotypes = matrix(nr=0,nc=0),
		groups = factor(0),
		sample_sizes = numeric(0),
		description = "",
		notes = ""
	)
)

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

