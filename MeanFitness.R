## CALCULATE MEAN AND VARIANCE FOR FITNESS

# for one population, within that pop, may later extend to do multiple pops


# NEED DIPLOID DATA
# so use full output at end


calc.fitness <- function(diploid.poly.muts.dat, full.genomes.dat, fixed.mut.dat, pop.size, generation){

	fitness.results <- data.frame(matrix(NA, ncol=2))
	names(fitness.results) <- c("individual", "poly.mut.fitness")	
	
	genome.data <- lapply(as.character(full.genomes.dat[,2]), FUN=strsplit, split=" ")

	# go through every two haploid genomes to get every diploid inds fitness
	iterate.inds <- 1
	for(i in seq(1, (pop.size*2), by=2)){
		genome.1.muts <- unlist(genome.data[[1]])[-1]
		genome.2.muts <- unlist(genome.data[[i+1]])[-1]
	
		# mutations both chromosomes have make them homozygotes for that mutation
		hom.muts <- genome.1.muts[which(genome.1.muts %in% genome.2.muts)]
		# mutations only one chromosome has make them heterozygotes for that mutation
		het.muts <- genome.1.muts[c(which(!(genome.1.muts %in% genome.2.muts)), which(!(genome.2.muts %in% genome.1.muts)))]
		# no need/no way to detect wild type homozygotes, so don't need to worry about for fitness
	
		
		## fitness is multiplicative:
		# remove neutral mutations because that will mess up the product
		hom.fitness <- prod(1 + diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% hom.muts) ,]$seln_coeff[diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% hom.muts) ,]$seln_coeff !=0])
		het.fitness <- prod(1 + (diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$seln_coeff[diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$seln_coeff !=0] * diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$dom_coeff[diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$seln_coeff !=0]))
	##	and take that value within this loop and save it, so that when getting decriment from fixed muts, can multiply together first then subtract from 1
				
		fitness.results[iterate.inds,] <- c(iterate.inds, (hom.fitness * het.fitness))
		iterate.inds <- iterate.inds + 1
	}
	
	# the above is all for polymorphic mutations
	
	# now need to include in individual fitness the effect of all fixed mutations
	#	(all fixed muts are always present at the last full generation time point sampled)
	fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]

	fixed.fitness <- prod(1 + fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff != 0])	# can't be right either because gives zero - so clearly all the realy neutral ones are the ones that are fixed, but they contribute nothing to fitness loss? Can't multiply them with the others because comes out to zero

	fitness.results$ind.fitness <- fixed.fitness * fitness.results$poly.mut.fitness
	
	
	mean.fitness.poly <- mean(fitness.results$poly.mut.fitness)
	var.fitness.poly <- var(fitness.results$poly.mut.fitness)
	mean.fitness.total <- mean(fitness.results$ind.fitness)
	var.fitness.total <- var(fitness.results$ind.fitness)
	
	return(c(mean.fitness.poly, var.fitness.poly, mean.fitness.total, var.fitness.total))
}



##	setwd("~/Documents/My_Documents/UofToronto/SLiM/AnalysisScripts")

##	generation <- 100000
##	num.inds.sampled <- 100
##	diploid.poly.muts.dat <- read.table("test_FULLfinal_polymuts.dat")
##	full.genomes.dat <- read.table("test_FULLfinal_genomes.dat", sep="A")
##	fixed.mut.dat <- read.table("test_fixedmuts_allgens.dat")

##	names(diploid.poly.muts.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
##	names(fixed.mut.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")


