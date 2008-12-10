if (!isGeneric("genotypes")) {
	if (is.function("genotypes"))
		fun <- genotypes
	else fun <- function(object) standardGeneric("genotypes")
	setGeneric("genotypes", fun)
}
	
if (!isGeneric("ploidy")) {
	if (is.function("ploidy"))
		fun <- ploidy
	else fun <- function(object) standardGeneric("ploidy")
	setGeneric("ploidy", fun)
}

if (!isGeneric("markers")) {
	if (is.function("markers"))
		fun <- markers
	else fun <- function(object) standardGeneric("markers")
	setGeneric("markers", fun)
}

if (!isGeneric("sampleSizes")) {
	if (is.function("sampleSizes"))
		fun <- sampleSizes
	else fun <- function(object) standardGeneric("sampleSizes")
	setGeneric("sampleSizes", fun)
}

if (!isGeneric("groups")) {
	if (is.function("groups"))
		fun <- groups
	else fun <- function(object) standardGeneric("groups")
	setGeneric("groups", fun)
}

if (!isGeneric("samples")) {
	if (is.function("samples"))
		fun <- samples
	else fun <- function(object) standardGeneric("samples")
	setGeneric("samples", fun)
}

if (!isGeneric("description")) {
	if (is.function("description"))
		fun <- description
	else fun <- function(object) standardGeneric("description")
	setGeneric("description", fun)
}

if (!isGeneric("notes")) {
	if (is.function("notes"))
		fun <- notes
	else fun <- function(object) standardGeneric("notes")
	setGeneric("notes", fun)
}

if (!isGeneric("extract")) {
	if (is.function("extract"))
		fun <- extract
	else fun <- function(object, ...) standardGeneric("extract")
	setGeneric("extract", fun)
}

