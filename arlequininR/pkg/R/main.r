"readArlequinData" <-
function (filename)
{	
	all <- scan(file=filename, what=character(), quiet=TRUE)

	## Find the Samples
	samples <- all[which(all == "[[Samples]]"):which(all == "[[Structure]]")]

	## xxx_i is xxx indices
	## xxx_e is xxx entries

	## Find the beginnings and ends of the samples
	## There will be an extra index in end_i for the closing bracket
	## 	for the file as a whole
	start_i <- which(samples == "{")
	end_i   <- which(samples == "}")

	## Find the sizes of each of the samples
	size_i  <- which(substr(samples, 1, 10) == "SampleSize")
	size_e  <- strsplit(samples[size_i], "=")
	sizes   <- c()
	for (i in 1:length(size_e)) {
		sizes <- c(sizes, as.numeric(size_e[[i]][2]))
	}

	## Within each sample, each individual has 2 entries for name and number
	##  of identical samples, followed by 2 entries for each locus listed as
	##  of the first entries then all of the second entries. 
	## So the number of loci is the number of slots for the sample divided
	##  by the number of individuals minus 2 divided by 2
	num_loci = ((end_i[1] - 1 - start_i[1])/sizes[1] - 2)/2

	## Pull out all the indices
	id_i        <- c()
	N_i         <- c()
	in_sample_i <- c()
	for (i in 1:length(sizes)) {
		id_i        <- c(id_i, 0:(sizes[i]-1)*(2+2*num_loci) + 1 + start_i[i])
		N_i         <- c(N_i,  0:(sizes[i]-1)*(2+2*num_loci) + 2 + start_i[i])
		in_sample_i <- c(in_sample_i, start_i[i]:end_i[i])
	}

	## Get rid of the ids, N values, and brackets for the data indices
	data_i <- setdiff(in_sample_i, c(id_i, N_i, start_i, end_i))

	## Pull out the sample ids, Ns, and data
	ids <- samples[id_i]
	Ns  <- as.numeric(samples[N_i])
	dat <- samples[data_i]

	sort_i <- c()
	for (i in 0:(sum(Ns)-1)) { 
		for (j in 1:num_loci) { 
			sort_i <- c(sort_i, i*2*num_loci + j, i*2*num_loci + num_loci + j)
		} 
	}

	## Use a generic marker identification scheme
	markers <- paste("M", rep(1:num_loci), sep="")

	## Put together the data matrix with named rows and columns
	dat_m <- matrix( dat[sort_i], 
					 nrow=sum(Ns), 
					 byrow=TRUE, 
					 dimnames=list(
						ids, 
						paste(rep(markers, each=2), "_", 1:2, sep="") ) )

	## Create the groups factor
	groups <- c()
	for (i in 1:length(sizes)) {
		groups <- c(groups, rep(i, sizes[i]))
	}
	groups <- factor(groups)

	s_sizes <- rep(1, length(ids))
	
	gd <- new( "genotypeData", 
			   ploidy=rep(2, num_loci), 
			   genotypes=dat_m, 
			   sample_sizes=s_sizes,
			   groups=groups,
			   markers=markers )
			  
	return(gd)
}

"conv" <-
function(a, b) 
{
  .C("convolve", as.double(a), as.integer(length(a)), as.double(b), as.integer(length(b)), ab = double(length(a) + length(b) - 1))$ab
}
