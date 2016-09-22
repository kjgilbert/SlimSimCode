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




calc.pi.stats <- function(mut.id.dat, genome.dat, num.inds.sampled, sequence.length){
	## WITH INFO ON NONSYNONYMOUS AND SYNONYMOUS MUTATIONS can also calculate pi_n/pi_s
	
	
	# make data frames of all possible neutral and deleterious mutations at all time points recorded
	neut.muts <- NULL
	seln.muts <- NULL
			
	## WHERE mut.id.dat IS A FILE OF ROWS OF MUTATIONS OCCURRING
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
	
	
	neut.muts <- mut.id.dat[mut.id.dat$mut.type == "m1" ,]
	seln.muts <- mut.id.dat[mut.id.dat$mut.type != "m1" ,]
			
	# all possible mut IDs of any type (need for later calcs):
	
	selected.mut.IDs <- c(seln.muts$mut.ID)
	neutral.mut.IDs <- c(neut.muts$mut.ID)
	
	
	## WHERE genome.dat IS A FILE OF MUTATIONS PER CHROMOSOME IN A SINGLE FILE LINE (IN SLIM, AFTER p1:532 A 1 2 3 4 5 ...)
	##		dat <- read.table(paste(c("genomes.out.", gens.sampled[gen], i), collapse=""), sep="A")
	##		column 1 = pop ID, then colon, then the ID of the individual, then A for Autosome and then the number of mutations with the identifiers 0 through ...
	
	
	
	# so take random pairs, find the number of non-overlapping mutations, i.e. unique, so polymorphic between them
	# sum across all pairs
	# divide by num genomes * length genome
	
	# get all possible line number pairs
	pairs <- combn(1:length(genome.dat[,1]), 2)
	total.poly <- 0
	nonsyn.total.poly <- 0
	syn.total.poly <- 0
	
	for(k in 1:dim(pairs)[2]){
		# take corresponding row for the pairs, and take second column of data file which has the mutation info
		pair1 <- genome.dat[pairs[1,k], 2]
		pair2 <- genome.dat[pairs[2,k], 2]
		
		total.muts1 <- unlist(strsplit(as.character(pair1), split=" "))[-1]	# must remove 1 b/c splitting the data makes a leading space that comes up as a shared mutation
		total.muts2 <- unlist(strsplit(as.character(pair2), split=" "))[-1]	# must remove 1 b/c splitting the data makes a leading space that comes up as a shared mutation
				
		neut.muts1 <- as.numeric(total.muts1[as.numeric(total.muts1) %in% neutral.mut.IDs])
		neut.muts2 <- as.numeric(total.muts2[as.numeric(total.muts2) %in% neutral.mut.IDs])
		seln.muts1 <- as.numeric(total.muts1[as.numeric(total.muts1) %in% selected.mut.IDs])
		seln.muts2 <- as.numeric(total.muts2[as.numeric(total.muts2) %in% selected.mut.IDs])
		
		# overall for pi
		total.matches <- length(intersect(total.muts1, total.muts2))	# how many mutations are the same
		poly1 <- length(total.muts1) - total.matches		# how many muts in ind 1 are not in ind 2
		poly2 <- length(total.muts2) - total.matches		# how many muts in ind 2 are not in ind 1
		temp.total.poly <- poly1 + poly2			# how many muts total are nonoverlapping between the two, i.e. the polymorphic/segregating sites
		total.poly <- total.poly + temp.total.poly
	
		# nonsynonymous (selected) for pi_n
		nonsyn.total.matches <- length(intersect(seln.muts1, seln.muts2))
		nonsyn.poly1 <- length(seln.muts1) - nonsyn.total.matches
		nonsyn.poly2 <- length(seln.muts2) - nonsyn.total.matches
		nonsyn.temp.total.poly <- nonsyn.poly1 + nonsyn.poly2
		nonsyn.total.poly <- nonsyn.total.poly + nonsyn.temp.total.poly
	
		# synonymous (neutral) for pi_s
		syn.total.matches <- length(intersect(neut.muts1, neut.muts2))
		syn.poly1 <- length(neut.muts1) - syn.total.matches
		syn.poly2 <- length(neut.muts2) - syn.total.matches
		syn.temp.total.poly <- syn.poly1 + syn.poly2
		syn.total.poly <- syn.total.poly + syn.temp.total.poly
	}
	
	pi <- total.poly / (num.inds.sampled * sequence.length)
	pi_n <- nonsyn.total.poly / (num.inds.sampled * sequence.length)
	pi_s <- syn.total.poly / (num.inds.sampled * sequence.length)

	return(c(pi, pi_n, pi_s))	
}