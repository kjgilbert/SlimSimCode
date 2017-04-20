
# segregating load - total (unrealized):

# segregating load - realized:

# fixed load:

calc.fixed.load <- function(fixeddat, from.gen, current.gen, pop.size){

	fixeddat <- fixeddat[fixeddat$gen.fixed > as.numeric(from.gen) & fixeddat$gen.fixed <= as.numeric(current.gen) ,]
		# this gives only mutations that have fixed PRIOR to and INCLUDING WITHIN the current generation time point sampled
	fixed.neut.muts <- c(which(fixeddat$mut.type == "m1"))
	fixed.seln.mut.IDs <- fixeddat$mut.ID[-fixed.neut.muts]
	fixed.neut.mut.IDs <- fixeddat$mut.ID[fixed.neut.muts]
			
	# break down selected fixed mutations into beneficial and deleterious, and Nes classes:
	seln.dat <- fixeddat[fixeddat$mut.ID %in% fixed.seln.mut.IDs ,]
	
	if(length(seln.dat$seln_coeff) == 0){
		total.fixed.load <- 0
	}else{
		total.fixed.load <- prod(seln.dat$seln_coeff, na.rm=TRUE)
	}
	if(length(seln.dat$seln_coeff[seln.dat$seln_coeff < 0]) == 0){
		total.del.load <- 0
	}else{
		total.del.load <- prod(seln.dat$seln_coeff[seln.dat$seln_coeff < 0], na.rm=TRUE)
	}
	if(length(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-1/pop.size) & seln.dat$seln_coeff < 0]) == 0){
		del.load_0_1 <- 0
	}else{
		del.load_0_1 <- prod(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-1/pop.size) & seln.dat$seln_coeff < 0], na.rm=TRUE)
	}
	if(length(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-10/pop.size) & seln.dat$seln_coeff < (-1/pop.size)]) == 0){
		del.load_1_10 <- 0
	}else{
		del.load_1_10 <- prod(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-10/pop.size) & seln.dat$seln_coeff < (-1/pop.size)], na.rm=TRUE)
	}
	if(length(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-100/pop.size) & seln.dat$seln_coeff < (-10/pop.size)]) == 0){
		del.load_10_100 <- 0
	}else{
		del.load_10_100 <- prod(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-100/pop.size) & seln.dat$seln_coeff < (-10/pop.size)], na.rm=TRUE)
	}
	if(length(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-Inf/pop.size) & seln.dat$seln_coeff < (-100/pop.size)]) == 0){
		del.load_100_inf <- 0
	}else{
		del.load_100_inf <- prod(seln.dat$seln_coeff[seln.dat$seln_coeff >= (-100/pop.size) & seln.dat$seln_coeff < (-10/pop.size)], na.rm=TRUE)
	}
	
	return(c(
		del.load_0_1,
		del.load_1_10,
		del.load_10_100,
		del.load_100_inf,
		total.del.load,
		total.fixed.load
		))
}

calc.seg.load <- function(poly.mut.dat, genome.data, from.gen, current.gen, num.inds){

	pop.total.load.w.recessive <- NULL
	pop.total.del.load.w.recessive <- NULL
	pop.total.load.w.recessive_0_1 <- NULL
	pop.total.load.w.recessive_1_10 <- NULL
	pop.total.load.w.recessive_10_100 <- NULL
	pop.total.load.w.recessive_100_inf <- NULL

	pop.total.load.w.realized <- NULL
	pop.total.del.load.w.realized <- NULL
	pop.total.load.w.realized_0_1 <- NULL
	pop.total.load.w.realized_1_10 <- NULL
	pop.total.load.w.realized_10_100 <- NULL
	pop.total.load.w.realized_100_inf <- NULL

	# go through inds, find muts per ind
	for(i in seq(1, (num.inds*2), by=2)){
		genome.1.muts <- as.numeric(unlist(strsplit(as.character(genome.data[[i,2]]), split=" "))[-1])
		genome.2.muts <- as.numeric(unlist(strsplit(as.character(genome.data[[(i+1),2]]), split=" "))[-1])
	
		# mutations both chromosomes have make them homozygotes for that mutation
		hom.muts <- genome.1.muts[which(genome.1.muts %in% genome.2.muts)]
		# mutations only one chromosome has make them heterozygotes for that mutation
		het.muts <- c(genome.1.muts[which(!(genome.1.muts %in% genome.2.muts))], genome.2.muts[which(!(genome.2.muts %in% genome.1.muts))])
		# no need/no way to detect wild type homozygotes, so don't need to worry about for fitness
		
		# take only the ones under seln, i.e. exclude s=0
		hom.muts.effect.sizes <- poly.mut.dat$seln_coeff[intersect(which(poly.mut.dat$mut.ID %in% hom.muts), which(poly.mut.dat$seln_coeff != 0))]
		het.muts.effect.sizes <- poly.mut.dat$seln_coeff[intersect(which(poly.mut.dat$mut.ID %in% het.muts), which(poly.mut.dat$seln_coeff != 0))]
		het.muts.dominance <- poly.mut.dat$dom_coeff[intersect(which(poly.mut.dat$mut.ID %in% het.muts), which(poly.mut.dat$seln_coeff != 0))]
		
	# RECESSIVE LOAD:		
		total.load.w.recessive.ind <- prod(hom.muts.effect.sizes) * prod(het.muts.effect.sizes)
		total.del.load.w.recessive.ind <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes < 0]) * prod(het.muts.effect.sizes[het.muts.effect.sizes < 0])
		total.load.w.recessive.ind_0_1 <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-1/pop.size) & hom.muts.effect.sizes < 0]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-1/pop.size) & het.muts.effect.sizes < 0])
		if(total.load.w.recessive.ind_0_1 == 1) total.load.w.recessive.ind_0_1 <- 0		# this would only happen if the subset was of length 0 then it took the product
		total.load.w.recessive.ind_1_10 <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-10/pop.size) & hom.muts.effect.sizes < (-1/pop.size)]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-10/pop.size) & het.muts.effect.sizes < (-1/pop.size)])
		if(total.load.w.recessive.ind_1_10 == 1) total.load.w.recessive.ind_1_10 <- 0		# this would only happen if the subset was of length 0 then it took the product
		total.load.w.recessive.ind_10_100 <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-100/pop.size) & hom.muts.effect.sizes < (-10/pop.size)]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-100/pop.size) & het.muts.effect.sizes < (-10/pop.size)])
		if(total.load.w.recessive.ind_10_100 == 1) total.load.w.recessive.ind_10_100 <- 0		# this would only happen if the subset was of length 0 then it took the product
		total.load.w.recessive.ind_100_inf <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-Inf/pop.size) & hom.muts.effect.sizes < (-100/pop.size)]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-Inf/pop.size) & het.muts.effect.sizes < (-100/pop.size)])
		if(total.load.w.recessive.ind_100_inf == 1) total.load.w.recessive.ind_100_inf <- 0		# this would only happen if the subset was of length 0 then it took the product
		
		pop.total.load.w.recessive <- c(pop.total.load.w.recessive, total.load.w.recessive.ind)
		pop.total.del.load.w.recessive <- c(pop.total.del.load.w.recessive, total.del.load.w.recessive.ind)
		pop.total.load.w.recessive_0_1 <- c(pop.total.load.w.recessive_0_1, total.load.w.recessive.ind_0_1)
		pop.total.load.w.recessive_1_10 <- c(pop.total.load.w.recessive_1_10, total.load.w.recessive.ind_1_10)
		pop.total.load.w.recessive_10_100 <- c(pop.total.load.w.recessive_10_100, total.load.w.recessive.ind_10_100)
		pop.total.load.w.recessive_100_inf <- c(pop.total.load.w.recessive_100_inf, total.load.w.recessive.ind_100_inf)

	# REALIZED LOAD:		
		total.load.w.realized.ind <- prod(hom.muts.effect.sizes) * prod(het.muts.effect.sizes*het.muts.dominance)
		total.del.load.w.realized.ind <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes < 0]) * prod(het.muts.effect.sizes[het.muts.effect.sizes < 0]*het.muts.dominance[het.muts.effect.sizes < 0])
		total.load.w.realized.ind_0_1 <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-1/pop.size) & hom.muts.effect.sizes < 0]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-1/pop.size) & het.muts.effect.sizes < 0]*het.muts.dominance[het.muts.effect.sizes >= (-1/pop.size) & het.muts.effect.sizes < 0])
		if(total.load.w.realized.ind_0_1 == 1) total.load.w.realized.ind_0_1 <- 0		# this would only happen if the subset was of length 0 then it took the product
		total.load.w.realized.ind_1_10 <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-10/pop.size) & hom.muts.effect.sizes < (-1/pop.size)]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-10/pop.size) & het.muts.effect.sizes < (-1/pop.size)]*het.muts.dominance[het.muts.effect.sizes >= (-10/pop.size) & het.muts.effect.sizes < (-1/pop.size)])
		if(total.load.w.realized.ind_1_10 == 1) total.load.w.realized.ind_1_10 <- 0		# this would only happen if the subset was of length 0 then it took the product
		total.load.w.realized.ind_10_100 <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-100/pop.size) & hom.muts.effect.sizes < (-10/pop.size)]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-100/pop.size) & het.muts.effect.sizes < (-10/pop.size)]*het.muts.dominance[het.muts.effect.sizes >= (-100/pop.size) & het.muts.effect.sizes < (-10/pop.size)])
		if(total.load.w.realized.ind_10_100 == 1) total.load.w.realized.ind_10_100 <- 0		# this would only happen if the subset was of length 0 then it took the product
		total.load.w.realized.ind_100_inf <- prod(hom.muts.effect.sizes[hom.muts.effect.sizes >= (-Inf/pop.size) & hom.muts.effect.sizes < (-100/pop.size)]) * prod(het.muts.effect.sizes[het.muts.effect.sizes >= (-Inf/pop.size) & het.muts.effect.sizes < (-100/pop.size)]*het.muts.dominance[het.muts.effect.sizes >= (-Inf/pop.size) & het.muts.effect.sizes < (-100/pop.size)])
		if(total.load.w.realized.ind_100_inf == 1) total.load.w.realized.ind_100_inf <- 0		# this would only happen if the subset was of length 0 then it took the product
	
		pop.total.load.w.realized <- c(pop.total.load.w.realized, total.load.w.realized.ind)
		pop.total.del.load.w.realized <- c(pop.total.del.load.w.realized, total.del.load.w.realized.ind)
		pop.total.load.w.realized_0_1 <- c(pop.total.load.w.realized_0_1, total.load.w.realized.ind_0_1)
		pop.total.load.w.realized_1_10 <- c(pop.total.load.w.realized_1_10, total.load.w.realized.ind_1_10)
		pop.total.load.w.realized_10_100 <- c(pop.total.load.w.realized_10_100, total.load.w.realized.ind_10_100)
		pop.total.load.w.realized_100_inf <- c(pop.total.load.w.realized_100_inf, total.load.w.realized.ind_100_inf)
	}
	total.load.w.recessive  <- mean(pop.total.load.w.recessive, na.rm=TRUE)
	total.del.load.w.recessive <- mean(pop.total.del.load.w.recessive, na.rm=TRUE)
	total.load.w.recessive_0_1 <- mean(pop.total.load.w.recessive_0_1, na.rm=TRUE)
	total.load.w.recessive_1_10 <- mean(pop.total.load.w.recessive_1_10, na.rm=TRUE)
	total.load.w.recessive_10_100 <- mean(pop.total.load.w.recessive_10_100, na.rm=TRUE)
	total.load.w.recessive_100_inf <- mean(pop.total.load.w.recessive_100_inf, na.rm=TRUE)
	
	total.load.w.realized <- mean(pop.total.load.w.realized, na.rm=TRUE)
	total.del.load.w.realized <- mean(pop.total.del.load.w.realized, na.rm=TRUE)
	total.load.w.realized_0_1 <- mean(pop.total.load.w.realized_0_1, na.rm=TRUE)
	total.load.w.realized_1_10 <- mean(pop.total.load.w.realized_1_10, na.rm=TRUE)
	total.load.w.realized_10_100 <- mean(pop.total.load.w.realized_10_100, na.rm=TRUE)
	total.load.w.realized_100_inf <- mean(pop.total.load.w.realized_100_inf, na.rm=TRUE)
	
	return(c(
	total.load.w.recessive,
	total.del.load.w.recessive,
	total.load.w.recessive_0_1,
	total.load.w.recessive_1_10,
	total.load.w.recessive_10_100,
	total.load.w.recessive_100_inf,
	total.load.w.realized,
	total.del.load.w.realized,
	total.load.w.realized_0_1,
	total.load.w.realized_1_10,
	total.load.w.realized_10_100,
	total.load.w.realized_100_inf	
	))
}





do.load.calcs <- function(sample.files, full.files=NULL, fixed.files, pop.size, inds.sampled, from.gen, current.gen, load.stats.output.file){

	results <- data.frame(matrix(nrow=0, ncol=21))
	names(results) <- c("ignore", "file", "generation", "del.load_0_1", "del.load_1_10", "del.load_10_100", "del.load_100_inf", "total.del.load", "total.fixed.load", "total.load.w.recessive", "total.del.load.w.recessive", "total.load.w.recessive_0_1", "total.load.w.recessive_1_10", "total.load.w.recessive_10_100", "total.load.w.recessive_100_inf", "total.load.w.realized", "total.del.load.w.realized", "total.load.w.realized_0_1", "total.load.w.realized_1_10", "total.load.w.realized_10_100", "total.load.w.realized_100_inf")
	write.table(results, append=FALSE, file=load.stats.output.file, sep=",", col.names=TRUE)

	for(i in 1:length(fixed.files)){
		# go through each file
		sample.file <- sample.files[i]
		fixed.file <- fixed.files[i]
		if(!is.null(full.files)){
			full.file <- full.files[i]
			## full data output
			full.samp.muts.start <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", full.file), collapse=""), intern=TRUE), split=":"))[1])
			full.samp.inds.start <- as.numeric(strsplit(system(paste(c("grep -n Individuals ", full.file), collapse=""), intern=TRUE), split=":")[[1]][1])
			full.samp.genomes.start <- as.numeric(strsplit(system(paste(c("grep -n Genomes ", full.file), collapse=""), intern=TRUE), split=":")[[1]][1])
			full.samp.file.end <- as.numeric(head(tail(unlist(strsplit(system(paste(c("wc -l ", full.file), collapse=""), intern=TRUE), split=" ")), n=2), n=1))
		}
	
		## sample data output
		poly.mut.id.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*(num.gens.sampled)), by=2)])
		generations <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n ' GS ' ", sample.file), collapse=""), intern=TRUE), split=" ")), ncol=4, byrow=TRUE)[,2])
		genomes.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Genomes ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*(num.gens.sampled)), by=2)])
	
		num.gens.sampled <- length(generations)
		if(!is.null(full.files)) num.gens.sampled <- num.gens.sampled + 1
		
		## fixed data output
		fixed.mut.id.start <- 2
		fixeddat <- read.table(fixed.file, skip=fixed.mut.id.start)
		names(fixeddat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")

		for(j in 1:num.gens.sampled){
			# go through each time point of a given file

			if(!is.null(full.files) & j == num.gens.sampled){
				gen <- format(j * pop.size, scientific=FALSE)
				
				polydat <- read.table(full.file, skip=full.samp.muts.start, nrow=((full.samp.inds.start-1) - full.samp.muts.start), sep=" ")
				names(polydat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")

				odd.nums <- seq(1, (pop.size * 2), by=2)
				sub.samp <- sample(odd.nums, size=inds.sampled, replace=FALSE)
				diploid.sub.samp <- sort(c(sub.samp, (sub.samp + 1)))

				genodat <- read.table(full.file, skip=full.samp.genomes.start, nrow=(pop.size*2), sep="A")
				genodat <- genodat[diploid.sub.samp ,]
			}else{
				gen <- generations[j]
				# for sample time points:
				polydat <- read.table(sample.file, skip=poly.mut.id.starts[j], nrow=((genomes.starts[j]-1) - poly.mut.id.starts[j]), sep=" ")
				names(polydat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
				genodat <- read.table(sample.file, skip=genomes.starts[j], nrow=(num.inds.sampled * 2), sep="A")		
			}
			fixed.load <- calc.fixed.load(fixeddat=fixeddat, from.gen, current.gen, pop.size=pop.size)
			seg.load <- calc.seg.load(poly.mut.dat=polydat, genome.data=genodat, from.gen, current.gen, num.inds=inds.sampled)
			
			temp.results <- c(sample.file, gen, fixed.load, seg.load)
			write.table(t(temp.results), append=TRUE, file=load.stats.output.file, sep=",", col.names=FALSE)
		}
		polydat <- NULL
		genodat <- NULL		
	}
}

