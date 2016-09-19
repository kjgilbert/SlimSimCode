## CALCULATE THE MEAN AND VARIANCE IN NUMBER OF DELETERIOUS MUTATIONS PER INDIVIDUAL



# do from each a file of sampled inds and their muts, then tack on data from fixed muts because that gives time mutation took to fix


mean.var.muts <- function(poly.mut.dat, ind.dat, generation, fixed.mut.dat, num.inds.sampled){
	
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
	delet.mut.IDs <- poly.mut.dat$mut.ID[-neut.muts]
	neut.mut.IDs <- poly.mut.dat$mut.ID[neut.muts]

	# then go through the genomes file and count how many of those mutations at that time point are in a given ind, do for all inds, then get mean and var
	genomes <- ind.dat[,2]
	
	muts.per.ind <- data.frame(matrix(NA, ncol=4))
	names(muts.per.ind) <- c("ind.ID", "num.delet.muts.poly", "num.neut.muts.poly", "total.muts.poly")

	for(i in 1:num.inds.sampled){
		num.delet.muts.ind <- table(as.integer(unlist(strsplit(as.character(genomes[i]), split=" "))) %in% delet.mut.IDs)[2]	# take only TRUEs
		num.neu.muts.ind <- table(as.integer(unlist(strsplit(as.character(genomes[i]), split=" "))) %in% neut.mut.IDs)[2]	# take only TRUEs
		muts.per.ind[i,] <- c(i, num.delet.muts.ind, num.neu.muts.ind, (num.delet.muts.ind + num.neu.muts.ind))
	}

	# the above is only for polymorphic mutations
	
	# now tack on numbers of fixed mutations	

	# only want at current generation or previous:
		
	fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= generation ,]
		# this gives only mutations that have fixed PRIOR to and INCLUDING WITHIN the current generation time point sampled
	fixed.neut.muts <- c(which(fixed.mut.dat$mut.type == "m1"))
	fixed.delet.mut.IDs <- fixed.mut.dat$mut.ID[-fixed.neut.muts]
	fixed.neut.mut.IDs <- fixed.mut.dat$mut.ID[fixed.neut.muts]
	
	num.neut.muts.fixed <- length(fixed.neut.mut.IDs)
	num.delet.muts.fixed <- length(fixed.delet.mut.IDs)

	# get the means and variances for polymorphic mutations of either type
	
	mean.delet.muts.per.ind.poly <- mean(muts.per.ind$num.delet.muts.poly)
	mean.neut.muts.per.ind.poly <- mean(muts.per.ind$num.neut.muts.poly)
	mean.total.muts.per.ind.poly <- mean(muts.per.ind$total.muts.poly)

	var.delet.muts.per.ind.poly <- var(muts.per.ind$num.delet.muts.poly)
	var.neut.muts.per.ind.poly <- var(muts.per.ind$num.neut.muts.poly)
	var.total.muts.per.ind.poly <- var(muts.per.ind$total.muts.poly)

	# include fixed mutations now:
	
	mean.delet.muts.per.ind.all <- mean((muts.per.ind$num.delet.muts.poly + num.delet.muts.fixed))
	mean.neut.muts.per.ind.all <- mean((muts.per.ind$num.neut.muts.poly + num.neut.muts.fixed))
	mean.total.muts.per.ind.all <- mean((muts.per.ind$total.muts.poly + num.delet.muts.fixed + num.neut.muts.fixed))

	var.delet.muts.per.ind.all <- var((muts.per.ind$num.delet.muts.poly + num.delet.muts.fixed))
	var.neut.muts.per.ind.all <- var((muts.per.ind$num.neut.muts.poly + num.neut.muts.fixed))
	var.total.muts.per.ind.all <- var((muts.per.ind$total.muts.poly + num.delet.muts.fixed + num.neut.muts.fixed))

	
	
	return(c(
		mean.delet.muts.per.ind.poly, var.delet.muts.per.ind.poly, 
		mean.neut.muts.per.ind.poly, var.neut.muts.per.ind.poly,
		mean.total.muts.per.ind.poly, var.total.muts.per.ind.poly,
		mean.delet.muts.per.ind.all, var.delet.muts.per.ind.all, 
		mean.neut.muts.per.ind.all, var.neut.muts.per.ind.all,
		mean.total.muts.per.ind.all, var.total.muts.per.ind.all,
		num.delet.muts.fixed, num.neut.muts.fixed))
}



##	setwd("~/Documents/My_Documents/UofToronto/SLiM/AnalysisScripts")

##	generation <- 20000
##	num.inds.sampled <- 100
##	poly.mut.dat <- read.table("test_polymuts_20000.dat")
##	ind.dat <- read.table("test_genomes_20000.dat", sep="A")
##	fixed.mut.dat <- read.table("test_fixedmuts_allgens.dat")

##	names(poly.mut.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
##	names(fixed.mut.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")







