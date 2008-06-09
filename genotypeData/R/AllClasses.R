setClass("genotypeData",
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
		ploidy = c(),
		markers = c(),
		genotypes = matrix(nr=0,nc=0),
		groups = c(),
		sample_sizes = c(),
		description = "",
		notes = ""
	)
)

