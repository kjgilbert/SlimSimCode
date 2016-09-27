## CALCULATE THE MEAN AND VARIANCE IN NUMBER OF DELETERIOUS MUTATIONS PER INDIVIDUAL



# do from each a file of sampled inds and their muts, then tack on data from fixed muts because that gives time mutation took to fix


mean.var.muts <- function(poly.mut.dat, genome.dat, generation, fixed.mut.dat, num.inds.sampled){
	
	# from sample data at a single generation and fixed data from one file for all generations
	
	# first go through poly.muts file at that time point and find the ID numbers for all the deleterious mutations

	# columns =
	#	(1) a unique identifying ID, 
	#	(2) the mutation type, 
	#	(3) the base position, 
	#	(4) the selection coefficient (here always 0 since this is a neutral model), 
	#	(5) the dominance coefficient (here always 0.5), 
	#	(6) the identifier of the subpopulation in which the mutation first arose, 
	#	(7) the generation in which it arose, and 
	#	(8) the prevalence of the mutation (the number of genomes that contain the mutation, where – in the way that SLiM uses the term “genome” – there are two genomes per individual)
	
	neut.muts <- c(which(poly.mut.dat$mut.type == "m1"))
	delet.muts <- c(which(poly.mut.dat$mut.type == "m2"), which(poly.mut.dat$mut.type == "m3"))
	ben.muts <- c(which(poly.mut.dat$mut.type == "m4"))
	delet.mut.IDs <- poly.mut.dat$mut.ID[delet.muts]
	ben.mut.IDs <- poly.mut.dat$mut.ID[ben.muts]
	neut.mut.IDs <- poly.mut.dat$mut.ID[neut.muts]

	# then go through the genomes file and count how many of those mutations at that time point are in a given ind, do for all inds, then get mean and var
	genomes <- genome.dat[,2]
	
	muts.per.ind <- data.frame(matrix(NA, ncol=5))
	names(muts.per.ind) <- c("ind.ID", "num.delet.muts.poly", "num.ben.muts.poly", "num.neut.muts.poly", "total.muts.poly")

	## THIS WAS FOR HAPLOIDS
	##	for(i in 1:num.inds.sampled){
	##		num.delet.muts.ind <- table(as.integer(unlist(strsplit(as.character(genomes[i]), split=" "))) %in% delet.mut.IDs)[2]	# take only TRUEs
	##		num.neu.muts.ind <- table(as.integer(unlist(strsplit(as.character(genomes[i]), split=" "))) %in% neut.mut.IDs)[2]	# take only TRUEs
	##		muts.per.ind[i,] <- c(i, num.delet.muts.ind, num.neu.muts.ind, (num.delet.muts.ind + num.neu.muts.ind))
	##	}
	
	## NOW FOR DIPLOIDS - need to get per IND, so over the two chromosomes
	iterate <- 1
	for(k in seq(1, (2* num.inds.sampled), by=2)){
		muts.chrom1 <- unlist(strsplit(as.character(genomes[k]), split=" "))[-1]
		muts.chrom2 <- unlist(strsplit(as.character(genomes[k+1]), split=" "))[-1]
		total.muts.ind <- unique(c(muts.chrom1, muts.chrom2))

		num.delet.muts.ind <- table(as.integer(total.muts.ind) %in% delet.mut.IDs)[2]	# take only TRUEs
		num.ben.muts.ind <- table(as.integer(total.muts.ind) %in% ben.mut.IDs)[2]	# take only TRUEs
		num.neut.muts.ind <- table(as.integer(total.muts.ind) %in% neut.mut.IDs)[2]		# take only TRUEs
		muts.per.ind[iterate,] <- c(iterate, num.delet.muts.ind, num.ben.muts.ind, num.neut.muts.ind, (num.delet.muts.ind + num.neut.muts.ind))
		iterate <- iterate + 1
	}


	# the above is only for polymorphic mutations
	
	# now tack on numbers of fixed mutations	

	# only want at current generation or previous:
		
	fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]
		# this gives only mutations that have fixed PRIOR to and INCLUDING WITHIN the current generation time point sampled
	fixed.neut.muts <- c(which(fixed.mut.dat$mut.type == "m1"))
	fixed.ben.muts <- c(which(fixed.mut.dat$mut.type == "m4"))
	fixed.delet.mut.IDs <- fixed.mut.dat$mut.ID[-c(fixed.neut.muts, fixed.ben.muts)]
	fixed.neut.mut.IDs <- fixed.mut.dat$mut.ID[fixed.neut.muts]
	fixed.ben.mut.IDs <- fixed.mut.dat$mut.ID[fixed.ben.muts]
	
	num.neut.muts.fixed <- length(fixed.neut.mut.IDs)
	num.ben.muts.fixed <- length(fixed.ben.mut.IDs)
	num.delet.muts.fixed <- length(fixed.delet.mut.IDs)

	# get the means and variances for polymorphic mutations of either type
	
	mean.delet.muts.per.ind.poly <- mean(muts.per.ind$num.delet.muts.poly)
	mean.ben.muts.per.ind.poly <- mean(muts.per.ind$num.ben.muts.poly)
	mean.neut.muts.per.ind.poly <- mean(muts.per.ind$num.neut.muts.poly)
	mean.total.muts.per.ind.poly <- mean(muts.per.ind$total.muts.poly)

	var.delet.muts.per.ind.poly <- var(muts.per.ind$num.delet.muts.poly)
	var.ben.muts.per.ind.poly <- var(muts.per.ind$num.ben.muts.poly)
	var.neut.muts.per.ind.poly <- var(muts.per.ind$num.neut.muts.poly)
	var.total.muts.per.ind.poly <- var(muts.per.ind$total.muts.poly)

	# include fixed mutations now:
	
	mean.delet.muts.per.ind.all <- mean((muts.per.ind$num.delet.muts.poly + num.delet.muts.fixed))
	mean.ben.muts.per.ind.all <- mean((muts.per.ind$num.ben.muts.poly + num.ben.muts.fixed))
	mean.neut.muts.per.ind.all <- mean((muts.per.ind$num.neut.muts.poly + num.neut.muts.fixed))
	mean.total.muts.per.ind.all <- mean((muts.per.ind$total.muts.poly + num.delet.muts.fixed + num.neut.muts.fixed))

	var.delet.muts.per.ind.all <- var((muts.per.ind$num.delet.muts.poly + num.delet.muts.fixed))
	var.ben.muts.per.ind.all <- var((muts.per.ind$num.ben.muts.poly + num.ben.muts.fixed))
	var.neut.muts.per.ind.all <- var((muts.per.ind$num.neut.muts.poly + num.neut.muts.fixed))
	var.total.muts.per.ind.all <- var((muts.per.ind$total.muts.poly + num.delet.muts.fixed + num.neut.muts.fixed))

	
	
	return(c(
		mean.delet.muts.per.ind.poly, var.delet.muts.per.ind.poly, 
		mean.ben.muts.per.ind.poly, var.ben.muts.per.ind.poly,
		mean.neut.muts.per.ind.poly, var.neut.muts.per.ind.poly,
		mean.total.muts.per.ind.poly, var.total.muts.per.ind.poly,
		mean.delet.muts.per.ind.all, var.delet.muts.per.ind.all, 
		mean.ben.muts.per.ind.all, var.ben.muts.per.ind.all, 
		mean.neut.muts.per.ind.all, var.neut.muts.per.ind.all,
		mean.total.muts.per.ind.all, var.total.muts.per.ind.all,
		num.delet.muts.fixed, num.ben.muts.fixed, num.neut.muts.fixed))
}



##	setwd("~/Documents/My_Documents/UofToronto/SLiM/AnalysisScripts")

##	generation <- 20000
##	num.inds.sampled <- 100
##	poly.mut.dat <- read.table("test_polymuts_20000.dat")
##	genome.dat <- read.table("test_genomes_20000.dat", sep="A")
##	fixed.mut.dat <- read.table("test_fixedmuts_allgens.dat")

##	names(poly.mut.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
##	names(fixed.mut.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")







