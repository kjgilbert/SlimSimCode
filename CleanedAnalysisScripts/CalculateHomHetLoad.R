

calc.hom.het.load.windows <- function(poly.mut.dat, full.genomes.dat, fixed.mut.dat, pop.size, generation, min.s=-1, max.s=0){
	# default window is all the negative s values

	genome.data <- lapply(as.character(full.genomes.dat[,2]), FUN=strsplit, split=" ")
	if(!is.null(fixed.mut.dat)) fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]


# for HETEROZYGOTES
		# get fitness in s windows
	het.window.del.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(het.window.del.fitnesses) <- c("individual", paste("s.window.het.delet.fitness", min.s, max.s, sep="_"))
		# get total fitness as well which includes beneficial mutations
	het.window.bendel.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(het.window.bendel.fitnesses) <- c("individual", paste("s.window.het.bendel.fitness", min.s, max.s, sep="_"))
		# also get numbers of mutations:
	het.window.num.del.muts <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(het.window.num.del.muts) <- c("individual", paste("s.window.het.delet.muts", min.s, max.s, sep="_"))
	het.window.num.bendel.muts <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(het.window.num.bendel.muts) <- c("individual", paste("s.window.het.bendel.muts", min.s, max.s, sep="_"))

# for HOMOZYGOTES
		# get fitness in s windows
	hom.window.del.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(hom.window.del.fitnesses) <- c("individual", paste("s.window.hom.delet.fitness", min.s, max.s, sep="_"))
		# get total fitness as well which includes beneficial mutations
	hom.window.bendel.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(hom.window.bendel.fitnesses) <- c("individual", paste("s.window.hom.bendel.fitness", min.s, max.s, sep="_"))
		# also get numbers of mutations:
	hom.window.num.del.muts <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(hom.window.num.del.muts) <- c("individual", paste("s.window.hom.delet.muts", min.s, max.s, sep="_"))
	hom.window.num.bendel.muts <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(hom.window.num.bendel.muts) <- c("individual", paste("s.window.hom.bendel.muts", min.s, max.s, sep="_"))

# for HOMOZYGOTES -- INCLUDE ALL FIXATIONS
		# get fitness in s windows
	homAndFixed.window.del.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(homAndFixed.window.del.fitnesses) <- c("individual", paste("s.window.homAndFixed.delet.fitness", min.s, max.s, sep="_"))
		# get total fitness as well which includes beneficial mutations
	homAndFixed.window.bendel.fitnesses <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(homAndFixed.window.bendel.fitnesses) <- c("individual", paste("s.window.homAndFixed.bendel.fitness", min.s, max.s, sep="_"))
		# also get numbers of mutations:
	homAndFixed.window.num.del.muts <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(homAndFixed.window.num.del.muts) <- c("individual", paste("s.window.homAndFixed.delet.muts", min.s, max.s, sep="_"))
	homAndFixed.window.num.bendel.muts <- data.frame(matrix(NA, ncol=(1+length(min.s))))
	names(homAndFixed.window.num.bendel.muts) <- c("individual", paste("s.window.homAndFixed.bendel.muts", min.s, max.s, sep="_"))
	
	
	iterate.inds <- 1
	total.bendel.num.het.muts <- 0
	total.bendel.num.hom.muts <- 0
	per.ind.bendel.num.het.muts <- NULL
	per.ind.bendel.num.hom.muts <- NULL
	total.del.num.het.muts <- 0
	total.del.num.hom.muts <- 0
	per.ind.del.num.het.muts <- NULL
	per.ind.del.num.hom.muts <- NULL
	
	for(k in seq(1, (pop.size*2), by=2)){			# go through individuals within generation within simulation replicate
		# go through every two haploid genomes to get every diploid inds fitness
		genome.1.muts <- unlist(genome.data[[k]])[-1]
		genome.2.muts <- unlist(genome.data[[k+1]])[-1]
	
		# mutations both chromosomes have make them homozygotes for that mutation
		hom.muts <- genome.1.muts[which(genome.1.muts %in% genome.2.muts)]
		# mutations only one chromosome has make them heterozygotes for that mutation
		het.muts <- c(genome.1.muts[which(!(genome.1.muts %in% genome.2.muts))], genome.2.muts[which(!(genome.2.muts %in% genome.1.muts))])
		# no need/no way to detect wild type homozygotes, so don't need to worry about for fitness
		
		# take only the ones under seln, i.e. exclude s=0
		hom.muts.effect.sizes <- poly.mut.dat$seln_coeff[intersect(which(poly.mut.dat$mut.ID %in% hom.muts), which(poly.mut.dat$seln_coeff != 0))]
		het.muts.effect.sizes <- poly.mut.dat$seln_coeff[intersect(which(poly.mut.dat$mut.ID %in% het.muts), which(poly.mut.dat$seln_coeff != 0))]
		het.muts.dominance <- poly.mut.dat$dom_coeff[intersect(which(poly.mut.dat$mut.ID %in% het.muts), which(poly.mut.dat$seln_coeff != 0))]



# at the end of k loop need total per pop and avg per ind:
#	num muts het and hom - total and per bin
#	fitness for het and hom muts - total and per bin
#	(2 categories of homs because can include or not include fixations in the pop from fixed data)

		window.num.del.het.muts <- NULL
		window.num.del.hom.muts <- NULL
		window.num.del.homAndFixed.muts <- NULL

		window.num.bendel.het.muts <- NULL
		window.num.bendel.hom.muts <- NULL
		window.num.bendel.homAndFixed.muts <- NULL
		
		window.het.del.fitness.per.ind <- NULL
		window.hom.del.fitness.per.ind <- NULL
		window.homAndFixed.del.fitness.per.ind <- NULL
		
		window.het.bendel.fitness <- NULL
		window.hom.bendel.fitness <- NULL
		window.homAndFixed.bendel.fitness <- NULL
				
		for(ll in 1:length(min.s)){			# go through mutations within individual within generation within simulation replicate to get fitness for those muts

			window.hom.fitness <- prod(1 + hom.muts.effect.sizes[hom.muts.effect.sizes > min.s[ll] & hom.muts.effect.sizes < max.s[ll]])	
			# because we're mostly doing negative s's, want the min to be listed as the largest negative number in the window
			window.het.fitness <- prod(1 + (het.muts.effect.sizes[het.muts.effect.sizes > min.s[ll] & het.muts.effect.sizes < max.s[ll]] * het.muts.dominance[het.muts.effect.sizes > min.s[ll] & het.muts.effect.sizes < max.s[ll]]))
			## fitness is multiplicative:
			window.het.del.fitness.per.ind <- c(window.het.del.fitness.per.ind, window.het.fitness)
			window.hom.del.fitness.per.ind <- c(window.hom.del.fitness.per.ind, window.hom.fitness)

			# with the beneficials:
			window.hom.fitness <- prod(1 + hom.muts.effect.sizes[c(intersect(which(hom.muts.effect.sizes > min.s[ll]), which(hom.muts.effect.sizes < max.s[ll])), which(hom.muts.effect.sizes > 0))])	
			# because we're mostly doing negative s's, want the min to be listed as the largest negative number in the window
			window.het.fitness <- prod(1 + (het.muts.effect.sizes[c(intersect(which(het.muts.effect.sizes > min.s[ll]), which(het.muts.effect.sizes < max.s[ll])), which(het.muts.effect.sizes > 0))] * het.muts.dominance[c(intersect(which(het.muts.effect.sizes > min.s[ll]), which(het.muts.effect.sizes < max.s[ll])), which(het.muts.effect.sizes > 0))]))
			## fitness is multiplicative:
			window.het.bendel.fitness <- c(window.het.bendel.fitness, window.het.fitness)
			window.hom.bendel.fitness <- c(window.hom.bendel.fitness, window.hom.fitness)

		# DO DELET MUTATION COUNTS on the current individual:
		
			window.num.del.hom.muts <- c(window.num.del.hom.muts, length(hom.muts.effect.sizes[hom.muts.effect.sizes > min.s[ll] & hom.muts.effect.sizes < max.s[ll]]))	
			# because we're mostly doing negative s's, want the min to be listed as the largest negative number in the window
			window.num.del.het.muts <- c(window.num.del.het.muts, length(het.muts.effect.sizes[het.muts.effect.sizes > min.s[ll] & het.muts.effect.sizes < max.s[ll]]))

			# with the beneficials:
			window.num.bendel.hom.muts <- c(window.num.bendel.hom.muts, length(hom.muts.effect.sizes[c(intersect(which(hom.muts.effect.sizes > min.s[ll]), which(hom.muts.effect.sizes < max.s[ll])), which(hom.muts.effect.sizes > 0))]))
			# because we're mostly doing negative s's, want the min to be listed as the largest negative number in the window
			window.num.bendel.het.muts <- c(window.num.bendel.het.muts, length(het.muts.effect.sizes[c(intersect(which(het.muts.effect.sizes > min.s[ll]), which(het.muts.effect.sizes < max.s[ll])), which(het.muts.effect.sizes > 0))]))
		}
		# get fitness for delet het muts in s windows
		het.window.del.fitnesses[iterate.inds ,] <- c(iterate.inds, window.het.del.fitness.per.ind)
		# get total fitness for het muts as well which includes beneficial mutations
		het.window.bendel.fitnesses[iterate.inds ,] <- c(iterate.inds, window.het.bendel.fitness)
		# also get numbers of heterozygous delet mutations:
		het.window.num.del.muts[iterate.inds ,] <- c(iterate.inds, window.num.del.het.muts)
		# also get numbers of heterozygous total mutations:
		het.window.num.bendel.muts[iterate.inds ,] <- c(iterate.inds, window.num.bendel.het.muts)

		# get fitness for delet hom muts in s windows
		hom.window.del.fitnesses[iterate.inds ,] <- c(iterate.inds, window.hom.del.fitness.per.ind)	
		# get total fitness for hom muts as well which includes beneficial mutations
		hom.window.bendel.fitnesses[iterate.inds ,] <- c(iterate.inds, window.hom.bendel.fitness)	
		# also get numbers of homozygous delet mutations:
		hom.window.num.del.muts[iterate.inds ,] <- c(iterate.inds, window.num.del.hom.muts)
		# also get numbers of homozygous total mutations:
		hom.window.num.bendel.muts[iterate.inds ,] <- c(iterate.inds, window.num.bendel.hom.muts)

		# also for HOMOZYGOTES - LATER - INCLUDE ALL FIXATIONS

		iterate.inds <- iterate.inds + 1
	}
	
	
	# the above is all for polymorphic mutations
	# now need to include in individual fitness the effect of all fixed mutations
	if(is.null(fixed.mut.dat)){			# but if nothing fixes:
		homAndFixed.window.del.fitnesses <- hom.window.del.fitnesses
		homAndFixed.window.bendel.fitnesses <- hom.window.bendel.fitnesses
		homAndFixed.window.num.del.muts <- hom.window.num.del.muts
		homAndFixed.window.num.bendel.muts <- hom.window.num.bendel.muts		
	}else{
		#	(all fixed muts are always present at the last full generation time point sampled)
		## already loaded it in at start of loop --- fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]
		
		window.fitnesses.fixed <- NULL
		window.fitnesses.wbens.fixed <- NULL
		# take numbers had and add new ones to these in the loop
		homAndFixed.window.num.del.muts <- hom.window.num.del.muts
		homAndFixed.window.num.bendel.muts <- hom.window.num.bendel.muts	
		for(ll in 1:length(min.s)){
			num.fixed.del.muts.in.bin <- length(fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > min.s[ll] & fixed.mut.dat$seln_coeff < max.s[ll]])
			num.fixed.bendel.muts.in.bin <- num.fixed.del.muts.in.bin + length(fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > 0])
			homAndFixed.window.num.del.muts[, ll+1] <- homAndFixed.window.num.del.muts[, ll+1] + num.fixed.del.muts.in.bin
			homAndFixed.window.num.bendel.muts[, ll+1] <- homAndFixed.window.num.bendel.muts[, ll+1] + num.fixed.bendel.muts.in.bin

			fixed.fitness.window <- prod(1 + fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > min.s[ll] & fixed.mut.dat$seln_coeff < max.s[ll]])
			window.fitnesses.fixed <- c(window.fitnesses.fixed, fixed.fitness.window)
			
			# totals with bens:
			fixed.fitness.wbens.window <- prod(1+fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > min.s[ll] & fixed.mut.dat$seln_coeff < max.s[ll] & fixed.mut.dat$seln_coeff > 0])
			window.fitnesses.wbens.fixed <- c(window.fitnesses.wbens.fixed, fixed.fitness.wbens.window)
		}
		names(window.fitnesses.fixed) <- paste("s.window.fixed.delet.fitness", min.s, max.s, sep="_")
		names(window.fitnesses.wbens.fixed) <- paste("s.window.fixed.bendel.fitness", min.s, max.s, sep="_")
		
		# total fitness is still multiplicative
		# combine with other homozygous fitnesses
		temp.hom <- hom.window.del.fitnesses[, -which(names(hom.window.del.fitnesses) == "individual")]
		homAndFixed.window.del.fitnesses <- data.frame(t(window.fitnesses.fixed * t(temp.hom)))

		temp.hom <- hom.window.bendel.fitnesses[, -which(names(hom.window.bendel.fitnesses) == "individual")]
		homAndFixed.window.bendel.fitnesses <- data.frame(t(window.fitnesses.wbens.fixed * t(temp.hom)))
	}

	mean.het.del.fitness <- apply(het.window.del.fitnesses[, - which(names(het.window.del.fitnesses) == "individual")], FUN=mean, MARGIN=2)
	mean.het.bendel.fitnesses <- apply(het.window.bendel.fitnesses[, - which(names(het.window.bendel.fitnesses) == "individual")], FUN=mean, MARGIN=2)
	mean.het.num.del.muts <- apply(het.window.num.del.muts[, - which(names(het.window.num.del.muts) == "individual")], FUN=mean, MARGIN=2)
	mean.het.num.bendel.muts <- apply(het.window.num.bendel.muts[, - which(names(het.window.num.bendel.muts) == "individual")], FUN=mean, MARGIN=2)

	mean.hom.del.fitness <- apply(hom.window.del.fitnesses[, - which(names(hom.window.del.fitnesses) == "individual")], FUN=mean, MARGIN=2)
	mean.hom.bendel.fitnesses <- apply(hom.window.bendel.fitnesses[, - which(names(hom.window.bendel.fitnesses) == "individual")], FUN=mean, MARGIN=2)
	mean.hom.num.del.muts <- apply(hom.window.num.del.muts[, - which(names(hom.window.num.del.muts) == "individual")], FUN=mean, MARGIN=2)
	mean.hom.num.bendel.muts <- apply(hom.window.num.bendel.muts[, - which(names(hom.window.num.bendel.muts) == "individual")], FUN=mean, MARGIN=2)

	mean.homAndFixed.del.fitness <- apply(homAndFixed.window.del.fitnesses[, - which(names(homAndFixed.window.del.fitnesses) == "individual")], FUN=mean, MARGIN=2)
	mean.homAndFixed.bendel.fitnesses <- apply(homAndFixed.window.bendel.fitnesses[, - which(names(homAndFixed.window.bendel.fitnesses) == "individual")], FUN=mean, MARGIN=2)
	mean.homAndFixed.num.del.muts <- apply(homAndFixed.window.num.del.muts[, - which(names(homAndFixed.window.num.del.muts) == "individual")], FUN=mean, MARGIN=2)
	mean.homAndFixed.num.bendel.muts <- apply(homAndFixed.window.num.bendel.muts[, - which(names(homAndFixed.window.num.bendel.muts) == "individual")], FUN=mean, MARGIN=2)
	
	return(c(
		mean.het.del.fitness, 
		mean.het.bendel.fitnesses,
		mean.het.num.del.muts, 
		mean.het.num.bendel.muts, 
		mean.hom.del.fitness, 
		mean.hom.bendel.fitnesses,
		mean.hom.num.del.muts, 
		mean.hom.num.bendel.muts, 
		mean.homAndFixed.del.fitness, 
		mean.homAndFixed.bendel.fitnesses,
		mean.homAndFixed.num.del.muts, 
		mean.homAndFixed.num.bendel.muts
	))
}
