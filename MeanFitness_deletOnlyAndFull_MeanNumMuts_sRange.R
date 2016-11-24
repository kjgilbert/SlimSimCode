## CALCULATE MEAN AND VARIANCE FOR FITNESS

# for one population, within that pop, may later extend to do multiple pops


# NEED DIPLOID DATA
# so use full output at end


calc.fitness.window <- function(poly.mut.dat, full.genomes.dat, fixed.mut.dat, pop.size, generation, min.s=-1, max.s=0){	# default window is all the negative s values

	genome.data <- lapply(as.character(full.genomes.dat[,2]), FUN=strsplit, split=" ")
	if(!is.null(fixed.mut.dat)) fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]
	
	# get fitness in s windows
	window.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(window.fitnesses) <- c("individual", paste("s.window.fitness.poly", min.s, max.s, sep="_"))

	# get total fitness as well which includes beneficial mutations
	window.fitnesses.wbens <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(window.fitnesses.wbens) <- c("individual", paste("s.window.fitnessAll.poly", min.s, max.s, sep="_"))
	
	iterate.inds <- 1
	for(i in seq(1, (pop.size*2), by=2)){
		# go through every two haploid genomes to get every diploid inds fitness
		genome.1.muts <- unlist(genome.data[[i]])[-1]
		genome.2.muts <- unlist(genome.data[[i+1]])[-1]
	
		# mutations both chromosomes have make them homozygotes for that mutation
		hom.muts <- genome.1.muts[which(genome.1.muts %in% genome.2.muts)]
		# mutations only one chromosome has make them heterozygotes for that mutation
		het.muts <- c(genome.1.muts[which(!(genome.1.muts %in% genome.2.muts))], genome.2.muts[which(!(genome.2.muts %in% genome.1.muts))])
		# no need/no way to detect wild type homozygotes, so don't need to worry about for fitness
		# length(unique(het.muts))
		
		# take only the ones under seln, i.e. exclude s=0
		hom.muts.effect.sizes <- poly.mut.dat[which(poly.mut.dat$mut.ID %in% hom.muts) ,]$seln_coeff[poly.mut.dat[which(poly.mut.dat$mut.ID %in% hom.muts) ,]$seln_coeff !=0]
		het.muts.effect.sizes <- poly.mut.dat[which(poly.mut.dat$mut.ID %in% het.muts) ,]$seln_coeff[poly.mut.dat[which(poly.mut.dat$mut.ID %in% het.muts) ,]$seln_coeff !=0]
		het.muts.dominance <- poly.mut.dat[which(poly.mut.dat$mut.ID %in% het.muts) ,]$dom_coeff[poly.mut.dat[which(poly.mut.dat$mut.ID %in% het.muts) ,]$seln_coeff !=0]
				
		window.fitnesses.per.ind <- NULL
		window.fitnessesAll.per.ind <- NULL	# for including total fitness, i.e. all beneficials at once on top of every window? - yes, don't break into positive s windows either
		for(j in 1:length(min.s)){
			window.hom.fitness <- prod(1 + hom.muts.effect.sizes[hom.muts.effect.sizes > min.s[j] & hom.muts.effect.sizes < max.s[j]])	
			# because we're mostly doing negative s's, want the min to be listed as the largest negative number in the window
			window.het.fitness <- prod(1 + (het.muts.effect.sizes[het.muts.effect.sizes > min.s[j] & het.muts.effect.sizes < max.s[j]] * het.muts.dominance[het.muts.effect.sizes > min.s[j] & het.muts.effect.sizes < max.s[j]]))
			## fitness is multiplicative:
			window.fitnesses.per.ind <- c(window.fitnesses.per.ind, (window.hom.fitness * window.het.fitness))

			# with the beneficials:
			window.hom.fitness <- prod(1 + hom.muts.effect.sizes[c(intersect(which(hom.muts.effect.sizes > min.s[j]), which(hom.muts.effect.sizes < max.s[j])), which(hom.muts.effect.sizes > 0))])	
			# because we're mostly doing negative s's, want the min to be listed as the largest negative number in the window
			window.het.fitness <- prod(1 + (het.muts.effect.sizes[c(intersect(which(het.muts.effect.sizes > min.s[j]), which(het.muts.effect.sizes < max.s[j])), which(het.muts.effect.sizes > 0))] * het.muts.dominance[c(intersect(which(het.muts.effect.sizes > min.s[j]), which(het.muts.effect.sizes < max.s[j])), which(het.muts.effect.sizes > 0))]))
			## fitness is multiplicative:
			window.fitnessesAll.per.ind <- c(window.fitnessesAll.per.ind, (window.hom.fitness * window.het.fitness))
			
			# DO DELET MUTATION COUNTS on the current individual:
			num.poly.delet.muts.in.this.window <- length( unique		poly.mut.dat$mut.ID[which(poly.mut.dat$mut.ID %in% c(genome.1.muts, genome.2.muts) & poly.mut.dat$seln_coeff]	)
			num.fixed.delet.muts.in.this.window
			num.total.delet.muts.in.this.window <- length(c(genome.1.muts, genome.2.muts)) # if I didn't need a window...
			
		}

		window.fitnesses[iterate.inds ,] <- c(iterate.inds, window.fitnesses.per.ind)
		##	and take that value within this loop and save it, so that when getting decriment from fixed muts, can multiply together first then subtract from 1
		window.fitnesses.wbens[iterate.inds ,] <- c(iterate.inds, window.fitnessesAll.per.ind)

		iterate.inds <- iterate.inds + 1
	}
	# the above is all for polymorphic mutations
		
	# now need to include in individual fitness the effect of all fixed mutations
		# but if nothing fixes:
	if(is.null(fixed.mut.dat)){
		fitness.results.polyANDfixed <- window.fitnesses		
		fitness.wbens.results.polyANDfixed <- window.fitnesses.wbens		
	}else{
		#	(all fixed muts are always present at the last full generation time point sampled)
		## already loaded it in at start of loop --- fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]
	
		window.fitnesses.fixed <- NULL
		window.fitnesses.wbens.fixed <- NULL
		for(j in 1:length(min.s)){
			fixed.fitness.window <- prod(1 + fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > min.s[j] & fixed.mut.dat$seln_coeff < max.s[j]])
			window.fitnesses.fixed <- c(window.fitnesses.fixed, fixed.fitness.window)
			# totals with bens:
			fixed.fitness.wbens.window <- prod(1 + fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > min.s[j] & fixed.mut.dat$seln_coeff < max.s[j] & fixed.mut.dat$seln_coeff > 0])
			window.fitnesses.wbens.fixed <- c(window.fitnesses.wbens.fixed, fixed.fitness.wbens.window)
		}
		# total fitness is still multiplicative
		poly.fitness <- window.fitnesses[, -which(names(window.fitnesses) == "individual")]
		fitness.results.polyANDfixed <- data.frame(t(window.fitnesses.fixed * t(poly.fitness)))
		names(fitness.results.polyANDfixed) <- paste("s.window.fitness.total", min.s, max.s, sep="_")

		# total fitness w/ bens is still multiplicative
		poly.fitness.wbens <- window.fitnesses.wbens[, -which(names(window.fitnesses.wbens) == "individual")]
		fitness.wbens.results.polyANDfixed <- data.frame(t(window.fitnesses.wbens.fixed * t(poly.fitness.wbens)))
		names(fitness.results.polyANDfixed) <- paste("s.window.fitness.wbens.total", min.s, max.s, sep="_")
	}
	
	mean.fitness.poly <- apply(window.fitnesses[, -which(names(window.fitnesses) == "individual")], FUN=mean, MARGIN=2)
	mean.fitness.total <- apply(fitness.results.polyANDfixed, FUN=mean, MARGIN=2)
	
	mean.fitness.wbens.poly <- apply(window.fitnesses.wbens[, -which(names(window.fitnesses.wbens) == "individual")], FUN=mean, MARGIN=2)
	mean.fitness.wbens.total <- apply(fitness.wbens.results.polyANDfixed, FUN=mean, MARGIN=2)

	
	return(c(mean.fitness.poly, mean.fitness.total, mean.fitness.wbens.poly, mean.fitness.wbens.total))
}


# returns number of windows*4

