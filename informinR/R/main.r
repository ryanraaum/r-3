"factor2numeric" <-
function (x)
{
	## Remove all entries with missing data
	y = na.omit(x)
	
	nrows <- dim(y)[1]
	ncols <- dim(y)[2]
	nmarkers <- y[1,"NumMarkers"]
	endinfo = which(names(y) == "NumMarkers") - 1
	startdata = endinfo + 2

	metadata <- y[seq(1,nrows,2),1:endinfo]	
	dat <- as.matrix(y[,startdata:ncols])
	
	## Make a transformation matrix to add successive rows
	tfm <- matrix(rep(diag(nrows/2), each=2), nrow=nrows/2, byrow=T)

	transformed <- tfm %*% dat
	
	## Remove all markers that do not vary
	transformed <- transformed[,setdiff(1:nmarkers, which(sd(transformed) == 0, arr.ind=T))]
	
	return(cbind(metadata, data.frame(NumMarkers=dim(transformed)[2]), transformed))
}

"pcaContribution" <-
function (x, prop=0.90)
{
	y <- factor2numeric(x)
	startdata <- which(names(y) == "NumMarkers") + 1
	nmarkers = y[1,"NumMarkers"]
	enddata   <- dim(y)[2]
	
	pr <- prcomp(y[,startdata:enddata], scale=T)
	proportions <- (pr$sdev)^2/sum((pr$sdev)^2)
	
	i <- 1
	while (sum(proportions[1:i]) < prop) {
		i <- i + 1
	}
	
	keep = pr$rotation[,1:i]
	
	score <- rep(0,nmarkers)
	
	for (j in 1:i) {
		multiplier <- 1/max(keep[,j])
		score <- score + abs(keep[,j]*multiplier)*(proportions[j]/proportions[1])
	}
	
	## Make the bottom end of the ranking 0.0
	score <- score - min(score)
	
	## Make the top end of the ranking 1.0
	score <- score / max(score)
	
	return(score)
}

