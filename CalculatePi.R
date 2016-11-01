## CALCULATE PI (nucleotide diversity)


#	nucleotide diversity, pi, in a sample of k homologous sequences (alleles)
#	the frequency at which randomly chosen pairs of alleles differ at a given nucleotide site
#	values for ind sites are averaged over all sites in the sequence to provide the summary stat for the whole sequence
# calculate via:
#	sum of numbers of differences between all possible pairs of alleles
#	divide by the number of sequence pairs that were compared and by the number of cases studied
#		(num seq pairs with k independent alleles, equals [k(k-1)]/2 )
#
# p23 C&C example on p 29
# 5022 bases long sequence
# 3 alleles (3 possible variants) so 3 pairwise comparisons
# 1 vs 2 differs at 7 sites, 1 vs 3 differs at 4 sites, 2 vs 3 differs at 3 sites
#	7+4+3 = 14
#	so pi = 14/(3*5022) = 0.00093




calc.pi.stats <- function(poly.dat, genome.dat, fixed.dat, generation, num.inds.sampled, use.manual.sample=FALSE){
	## WITH INFO ON NONSYNONYMOUS AND SYNONYMOUS MUTATIONS can also calculate pi_n/pi_s
	
	# because diploid:
	sample.size <- 2 * num.inds.sampled
	
	# make data frames of all possible neutral and deleterious mutations at all time points recorded
	neut.muts <- NULL
	seln.muts <- NULL

	# tack on fixed data for the counts of fixed alleles not just in the sample but across the whole pop
	if(is.null(fixed.dat)){	# unless nothing has fixed
		num.neut.muts.fixed <- 0
		num.seln.muts.fixed <- 0
	}else{
		fixed.mut.dat <- fixed.dat[fixed.dat$gen.fixed <= as.numeric(generation) ,]
			# this gives only mutations that have fixed PRIOR to and INCLUDING WITHIN the current generation time point sampled
		fixed.neut.muts <- c(which(fixed.mut.dat$mut.type == "m1"))
		fixed.seln.mut.IDs <- fixed.mut.dat$mut.ID[-fixed.neut.muts]
		fixed.neut.mut.IDs <- fixed.mut.dat$mut.ID[fixed.neut.muts]
		
		num.neut.muts.fixed <- length(fixed.neut.mut.IDs)
		num.seln.muts.fixed <- length(fixed.seln.mut.IDs)
	}
			
	## WHERE poly.dat IS A FILE OF ROWS OF MUTATIONS OCCURRING
	#	m1 = neutral site in coding
	#	m2, 3 = selected site in coding
	
	# columns =
	#	(1) a unique identifying ID, 
	#	(2) the mutation type, 
	#	(3) the base position, 
	#	(4) the selection coefficient (here always 0 since this is a neutral model), 
	#	(5) the dominance coefficient (here always 0.5), 
	#	(6) the identifier of the subpopulation in which the mutation first arose, 
	#	(7) the generation in which it arose, and 
	#	(8) the prevalence of the mutation (the number of genomes that contain the mutation, where – in the way that SLiM uses the term “genome” – there are two genomes per individual)
	
	##	dat <- read.table(paste(c("poly.muts.out.", gens.sampled[gen], i), collapse=""))
	##	names(dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
	
	
	neut.muts <- poly.dat[poly.dat$mut.type == "m1" ,]
	seln.muts <- poly.dat[poly.dat$mut.type != "m1" ,]
			
	if(use.manual.sample == FALSE){
		sfs.total <- table(poly.dat$mut.prev)
		sfs.neut <- table(neut.muts$mut.prev)
		sfs.seln <- table(seln.muts$mut.prev)
	}
	if(use.manual.sample == TRUE){
		## would have to use this section of code when manually subsampling a full sample output, because don't have the allele frequencies in the sample...
		
		# all possible mut IDs of any type (need for later calcs):
		
		neutral.mut.IDs <- c(neut.muts$mut.ID)
		selected.mut.IDs <- c(seln.muts$mut.ID)
		
		
		all.mutations <- unlist(lapply(as.character(genome.dat[,2]), FUN=strsplit, split=" "))
		just.neut.muts <- all.mutations[which(all.mutations %in% neutral.mut.IDs)]
		just.seln.muts <- all.mutations[which(all.mutations %in% selected.mut.IDs)]
		
		# get the allele frequencies:
		freqs.total <- table(all.mutations)
		freqs.neut <- table(just.neut.muts)
		freqs.seln <- table(just.seln.muts)

		# get the frequency spectra
		sfs.total <- table(freqs.total)
		sfs.neut <- table(freqs.neut)
		sfs.seln <- table(freqs.seln)
	}
	
## NO don't add fixed stuff because it's supposed to measure POLYMORPHISM and polymorphic sites! does end up mattering b/c I divide by all counts (the 2pq calcs go to zero though)
##	# add on fixed things to the fixed section of the table
##	if(is.na(sfs.total[as.character(sample.size)])){
##		sfs.total[as.character(sample.size)] <- num.neut.muts.fixed + num.seln.muts.fixed
##	}else{
##		sfs.total[as.character(sample.size)] <- sfs.total[as.character(sample.size)] + num.neut.muts.fixed + num.seln.muts.fixed
##	}
##	if(is.na(sfs.neut[as.character(sample.size)])){
##		sfs.neut[as.character(sample.size)] <- num.neut.muts.fixed
##	}else{
##		sfs.neut[as.character(sample.size)] <- sfs.total[as.character(sample.size)] + num.neut.muts.fixed
##	}
##	if(is.na(sfs.seln[as.character(sample.size)])){
##		sfs.seln[as.character(sample.size)] <- num.seln.muts.fixed
##	}else{
##		sfs.seln[as.character(sample.size)] <- sfs.seln[as.character(sample.size)] + num.seln.muts.fixed
##	}


	counts.sfs.total <- as.numeric(names(sfs.total))
	p.all <- counts.sfs.total/sample.size
	q.all <- 1-p.all
	numerator.all <- 2*p.all*q.all*sfs.total
	pi.all <- sum(numerator.all)/sum(sfs.total)

	counts.sfs.neut <- as.numeric(names(sfs.neut))
	p.neut <- counts.sfs.neut/sample.size
	q.neut <- 1-p.neut
	numerator.neut <- 2*p.neut*q.neut*sfs.neut
	pi.synonymous <- sum(numerator.neut)/sum(sfs.neut)

	counts.sfs.seln <- as.numeric(names(sfs.seln))
	p.seln <- counts.sfs.seln/sample.size
	q.seln <- 1-p.seln
	numerator.seln <- 2*p.seln*q.seln*sfs.seln
	pi.nonsynonymous <- sum(numerator.seln)/sum(sfs.seln)
	
	pi <- pi.all
	pi_n <- pi.nonsynonymous
	pi_s <- pi.synonymous

	return(c(pi, pi_n, pi_s))	
}