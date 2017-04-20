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
# so theta_w = 7 / (5022*1.5) = 0.00093			(same as pi, which is expected if pop is at equil and no sites are under seln)






# harmonic sum function for theta
sum.a_k <- function(k){
	harmsum <- 0 
	for(i in 1:(k-1)){
		harmsum <- harmsum + 1/i
		if (harmsum >= (k-1)) break
	}
	return(harmsum)
}


## WHERE genome.dat IS A FILE OF MUTATIONS PER CHROMOSOME IN A SINGLE LINE (IN SLIM, AFTER p1:532 A 1 2 3 4 5 ...)
##		genome.dat <- read.table(paste(c("genomes.out.", gens.sampled[gen], i), collapse=""), sep="A")
##		column 1 = pop ID, then colon, then the ID of the individual, then A for Autosome and then the number of mutations with the identifiers 0 through ...

calc.theta <- function(genome.dat, poly.dat, generation, num.inds.sampled, sequence.length){
	# get total num poly sites across all samples = S_k; i.e., how many muts are not fixed in the sample?
	
	# because it's diploids, sampling 100 inds gives 200 sequences, so:
	sample.size <- 2 * num.inds.sampled
	
	# first pair shared mutations:
	geno1muts <- unlist(strsplit(as.character(genome.dat[1,2]), split=" "))
	geno2muts <- unlist(strsplit(as.character(genome.dat[2,2]), split=" "))
	matched.sites <- intersect(geno1muts, geno2muts)
	
	# get neutral muts only and just compare those:
	neut.muts.poly <- poly.dat$mut.ID[poly.dat$mut.type == "m1"]
	matched.sites.neut <- intersect(matched.sites, neut.muts.poly)
	
	# those shared mutations with 3rd individual and onward:
	for(k in 3:sample.size){
		next.geno.muts <- unlist(strsplit(as.character(genome.dat[k,2]), split=" "))
		temp.matched.sites <- intersect(next.geno.muts, matched.sites)

		matched.sites <- temp.matched.sites
	}
	matched.sites.neut <- matched.sites[matched.sites %in% as.character(neut.muts.poly)]
	
	num.matched.sites <- length(matched.sites) - 1	# must subtract 1 because when I split the genome.data into the second part, there is a leading space, so when it's split again, the space comes up as a shared mutation
	num.matched.sites.neut <- length(matched.sites.neut)	# no subtraction because this is generated in a way to remove the empty slot in the vector
	
	# that many ARE fixed in the SAMPLE, so total number muts in sample minus this should be number polymorphic sites in the sample
	
	# to get S_k, the total number of polymorphic sites, count all mutations and subtract the number matched across all
	all.muts <- 0
	for(k in 1:sample.size){
		temp.all.muts <- unlist(strsplit(as.character(genome.dat[k,2]), split=" "))
		all.muts <- c(all.muts, temp.all.muts)
	}
	all.muts.neut <- all.muts[all.muts %in% as.character(neut.muts.poly)]
	
	num.all.muts <- length(unique(all.muts)) - 1
		# must remove 1 b/c splitting the genome.data makes a leading space that comes up as a shared mutation
	num.all.muts.neut <- length(unique(all.muts.neut))	# no - 1 b/c generated in a way to remove the empty slot in the vector
	
	num.poly.muts <- num.all.muts - num.matched.sites
	num.poly.muts.neut <- num.all.muts.neut - num.matched.sites.neut
		
	theta <- num.poly.muts / (sum.a_k(sample.size) * sequence.length)
	theta.neut <- num.poly.muts.neut / (sum.a_k(sample.size) * (sequence.length*0.25))
	
	return(c(theta, theta.neut))
}