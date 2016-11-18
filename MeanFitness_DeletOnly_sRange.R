## CALCULATE MEAN AND VARIANCE FOR FITNESS

# for one population, within that pop, may later extend to do multiple pops


# NEED DIPLOID DATA
# so use full output at end


calc.fitness.window <- function(diploid.poly.muts.dat, full.genomes.dat, fixed.mut.dat, pop.size, generation, min.s=-1, max.s=0){	# default window is all the negative s values

	genome.data <- lapply(as.character(full.genomes.dat[,2]), FUN=strsplit, split=" ")

	# get fitness in s windows
	window.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(window.fitnesses) <- c("individual", paste("s.window", min.s, max.s, sep="_"))
	
	iterate.inds <- 1
	for(i in seq(1, (pop.size*2), by=2)){
		# go through every two haploid genomes to get every diploid inds fitness
		genome.1.muts <- unlist(genome.data[[i]])[-1]
		genome.2.muts <- unlist(genome.data[[i+1]])[-1]
	
		# mutations both chromosomes have make them homozygotes for that mutation
		hom.muts <- genome.1.muts[which(genome.1.muts %in% genome.2.muts)]
		# mutations only one chromosome has make them heterozygotes for that mutation
		het.muts <- genome.1.muts[c(which(!(genome.1.muts %in% genome.2.muts)), which(!(genome.2.muts %in% genome.1.muts)))]
		# no need/no way to detect wild type homozygotes, so don't need to worry about for fitness
		# length(unique(het.muts))
		
		# take only the ones under seln, i.e. exclude s=0
		hom.muts.effect.sizes <- diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% hom.muts) ,]$seln_coeff[diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% hom.muts) ,]$seln_coeff !=0]
		het.muts.effect.sizes <- diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$seln_coeff[diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$seln_coeff !=0]
		het.muts.dominance <- diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$dom_coeff[diploid.poly.muts.dat[which(diploid.poly.muts.dat$mut.ID %in% het.muts) ,]$seln_coeff !=0]
				
		window.fitnesses.per.ind <- NULL
		for(j in 1:length(min.s)){
			window.hom.fitness <- prod(1 + hom.muts.effect.sizes[hom.muts.effect.sizes > min.s[j] & hom.muts.effect.sizes < max.s[j]])	
			# because we're mostly doing negative s's, want the min to be listed as the largest negative number in the window
			window.het.fitness <- prod(1 + (het.muts.effect.sizes[het.muts.effect.sizes > min.s[j] & het.muts.effect.sizes < max.s[j]] * het.muts.dominance[het.muts.effect.sizes > min.s[j] & het.muts.effect.sizes < max.s[j]]))
			## fitness is multiplicative:
			window.fitnesses.per.ind <- c(window.fitnesses.per.ind, (window.hom.fitness * window.het.fitness))
		}
		window.fitnesses[iterate.inds ,] <- c(iterate.inds, window.fitnesses.per.ind)
		##	and take that value within this loop and save it, so that when getting decriment from fixed muts, can multiply together first then subtract from 1
		iterate.inds <- iterate.inds + 1
	}
	# the above is all for polymorphic mutations
	
	# now need to include in individual fitness the effect of all fixed mutations
		# but if nothing fixes:
	if(is.null(fixed.mut.dat)){
		fitness.results$ind.fitness <- fitness.results$poly.mut.fitness		
	}else{
		#	(all fixed muts are always present at the last full generation time point sampled)
		fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]
	
		window.fitnesses.fixed <- NULL
		for(j in 1:length(min.s)){
			fixed.fitness.window <- prod(1 + fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > min.s[j] & fixed.mut.dat$seln_coeff < max.s[j]])
			window.fitnesses.fixed <- c(window.fitnesses.fixed, fixed.fitness.window)
		}
		# total fitness is still multiplicative
		fitness.results.polyANDfixed <- window.fitnesses.fixed * window.fitnesses
	}
	
	mean.fitness.poly <- colMeans(window.fitnesses[, -which(names(window.fitnesses) == "individual")])
	var.fitness.poly <- apply(window.fitnesses[, -which(names(window.fitnesses) == "individual")], FUN=var, MARGIN=2)
	mean.fitness.total <- colMeans(fitness.results.polyANDfixed[, -which(names(fitness.results.polyANDfixed) == "individual")])
	var.fitness.total <- apply(fitness.results.polyANDfixed[, -which(names(fitness.results.polyANDfixed) == "individual")], FUN=var, MARGIN=2)
	
	return(c(mean.fitness.poly, var.fitness.poly, mean.fitness.total, var.fitness.total))
}


# returns number of windows*4

