## CALCULATE WATTERSON'S THETA





#	Watterson's theta
#	k = number of sequences in the sample
#	counts number of polymorphic sites in the entire sample = S_k
#	divide by a_k, which = harmonic mean over k
#		if k is large, a_k ~= ln(k)
#	theta_w = S_k / (a_k * number of bases in sequence)
#	a_k = 1 + 1/2 + 1/3 + 1/4 + ... + 1/(k-1)
# k = 3, so a_k = 1.5
# S_k = 7 because 7 total sites differ mong any of the 3 sequences
# length = 5022
# so theta_w = 7 / (5022*1.5) = 0.00093			(same as pi, which is expected if pop is at equil)






# harmonic sum function for theta
sum.a_k <- function(k){
	harmsum <- 0 
	for(i in 1:(k-1)){
		harmsum <- harmsum + 1/i
		if (harmsum >= (k-1)) break
	}
	return(harmsum)
}


## WHERE DAT IS A FILE OF MUTATIONS PER CHROMOSOME IN A SINGLE FILE LINE (IN SLIM, AFTER p1:532 A 1 2 3 4 5 ...)
##		dat <- read.table(paste(c("genomes.out.", gens.sampled[gen], i), collapse=""), sep="A")
##		column 1 = pop ID, then colon, then the ID of the individual, then A for Autosome and then the number of mutations with the identifiers 0 through ...


# get total num poly sites across all samples = S_k; i.e., how many muts are not fixed in the sample?

# first pair shared mutations:
matched.sites <- intersect(unlist(strsplit(as.character(dat[1,2]), split=" ")), unlist(strsplit(as.character(dat[2,2]), split=" ")))

# those shared mutations with 3rd individual and onward:
for(k in 3:num.inds.sampled){
	temp.matched.sites <- intersect(unlist(strsplit(as.character(dat[k,2]), split=" ")), matched.sites)
	matched.sites <- temp.matched.sites
}

num.matched.sites <- length(matched.sites) - 1	# must subtract 1 because when I split the data into the second part, there is a leading space, so when it's split again, the space comes up as a shared mutation

# that many ARE fixed in the SAMPLE, so total number muts in sample minus this should be number polymorphic sites in the sample

# to get S_k, the total number of polymorphic sites, count all mutations and subtract the number matched across all
all.muts <- 0
for(k in 1:num.inds.sampled){
	temp.all.muts <- unlist(strsplit(as.character(dat[k,2]), split=" "))
	all.muts <- c(all.muts, temp.all.muts)
}

num.all.muts <- length(unique(all.muts)) - 1
	# must remove 1 b/c splitting the data makes a leading space that comes up as a shared mutation

num.poly.muts <- num.all.muts - num.matched.sites
	
theta <- num.poly.muts / (sum.a_k(num.inds.sampled) * sequence.length)
