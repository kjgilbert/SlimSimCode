## CALCULATE DN/DS per effect size class of mutations




dnds.summ.stats <- function(fixed.output.files, outfile.name, start.gens, end.gens, geno.size, pop.size){
	
	results <- data.frame(matrix(nrow=0, ncol=11))
	names(results) <- c("ignore", "file", "generation.start", "generation.stop", "dnds_nes_0_1", "dnds_nes_1_10", "dnds_nes_10_100", "dnds_nes_100_inf", "dnds_totaldelet", "dnds_totalben", "dnds_totaloverall")
	write.table(results, append=FALSE, file=outfile.name, sep=",", col.names=TRUE)
	
	# gen start and stop is the span of generations over which dnds is calculated in the file

	# this is the subfunction for calculating dnds
	calc.dn.ds <- function(fixeddat, geno.size, from.gen, current.gen){
	
		fixeddat <- fixeddat[fixeddat$gen.fixed > as.numeric(from.gen) & fixeddat$gen.fixed <= as.numeric(current.gen) ,]
			# this gives only mutations that have fixed PRIOR to and INCLUDING WITHIN the current generation time point sampled
		fixed.neut.muts <- c(which(fixeddat$mut.type == "m1"))
		fixed.seln.mut.IDs <- fixeddat$mut.ID[-fixed.neut.muts]
		fixed.neut.mut.IDs <- fixeddat$mut.ID[fixed.neut.muts]
		
		num.neut.muts.fixed <- length(fixed.neut.mut.IDs)
		num.seln.muts.fixed <- length(fixed.seln.mut.IDs)
		
		# break down selected fixed mutations into beneficial and deleterious, and Nes classes:
		seln.dat <- fixeddat[fixeddat$mut.ID %in% fixed.seln.mut.IDs ,]
		num.ben.muts <- table(seln.dat$seln_coeff > 0)["TRUE"]
		if(is.na(num.ben.muts)) num.ben.muts <- 0
	
		num.del.muts <- table(seln.dat$seln_coeff < 0)["TRUE"]
		if(is.na(num.del.muts)) num.del.muts <- 0
	
		num.del.muts_0_1 <- table(seln.dat$seln_coeff >= (-1/pop.size) & seln.dat$seln_coeff < 0)["TRUE"]
		if(is.na(num.del.muts_0_1)) num.del.muts_0_1 <- 0
		num.del.muts_1_10 <- table(seln.dat$seln_coeff >= (-10/pop.size) & seln.dat$seln_coeff < (-1/pop.size))["TRUE"]
		if(is.na(num.del.muts_1_10)) num.del.muts_1_10 <- 0
		num.del.muts_10_100 <- table(seln.dat$seln_coeff >= (-100/pop.size) & seln.dat$seln_coeff < (-10/pop.size))["TRUE"]
		if(is.na(num.del.muts_10_100)) num.del.muts_10_100 <- 0
		num.del.muts_100_inf <- table(seln.dat$seln_coeff >= (-Inf/pop.size) & seln.dat$seln_coeff < (-100/pop.size))["TRUE"]
		if(is.na(num.del.muts_100_inf)) num.del.muts_100_inf <- 0
	
		# dn ds = selected substitutions per genome size / neutral substitutions per genome size
		dnds_nes_0_1 <- (num.del.muts_0_1/(geno.size*0.75)) / (num.neut.muts.fixed/(geno.size*0.25))
		dnds_nes_1_10 <- (num.del.muts_1_10/(geno.size*0.75)) / (num.neut.muts.fixed/(geno.size*0.25))
		dnds_nes_10_100 <- (num.del.muts_10_100/(geno.size*0.75)) / (num.neut.muts.fixed/(geno.size*0.25))
		dnds_nes_100_inf <- (num.del.muts_100_inf/(geno.size*0.75)) / (num.neut.muts.fixed/(geno.size*0.25))
		dnds_totaldelet <- (num.del.muts/(geno.size*0.75)) / (num.neut.muts.fixed/(geno.size*0.25))
		dnds_totalben <- (num.ben.muts/(geno.size*0.75)) / (num.neut.muts.fixed/(geno.size*0.25))
		dnds_totaloverall <- ((num.del.muts + num.ben.muts)/(geno.size*0.75)) / (num.neut.muts.fixed/(geno.size*0.25))
		
		return(c(
			dnds_nes_0_1,
			dnds_nes_1_10,
			dnds_nes_10_100,
			dnds_nes_100_inf,
			dnds_totaldelet,
			dnds_totalben,
			dnds_totaloverall
			))
	}
			
	iterate <- 1
	for(i in 1:length(fixed.output.files)){
		# go through each file
		
		fixed.file <- fixed.output.files[i]
		
		## fixed data output
		if(length(readLines(fixed.file)) == 2){	# then no mutations fixed
			fixeddat <- NULL
		}else{	# otherwise read in fixed mutations as normal
			fixed.mut.id.start <- 2
			fixeddat <- read.table(fixed.file, skip=fixed.mut.id.start)
			names(fixeddat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
		}
	
		for(j in 1:length(start.gens)){
			from.gen <- start.gens[j]
			current.gen <- end.gens[j]
			
			dnds.stats  <- calc.dn.ds(fixeddat=fixeddat, geno.size=geno.size, from.gen=from.gen, current.gen=current.gen)
					
			temp.results <- c(fixed.file, from.gen, current.gen, dnds.stats)
			write.table(t(temp.results), append=TRUE, file=outfile.name, sep=",", col.names=FALSE)			
		}
		
		iterate <- iterate + 1
		fixeddat <- NULL
	}
	print("Calculations complete :) ")
}





