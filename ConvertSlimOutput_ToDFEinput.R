
# convert last time point output to SFS:
#
#	The folded frequency spectrum stores the observed counts of the minor (most rare) allele frequencies. The folded spectrum can be calculated by binning together the ith and (n-i)th entries from the unfolded spectrum, where n is the number of sampled chromosomes.
#	The program to estimate the DFE is called est_dfe, and can be run in two modes, depending on whether the folded or unfolded SFS is analysed. It is more straightforward and the results are likely to be more robust from analysis of the folded SFS.
#	If the folded spectra are analyzed, running the program est_alpha_omega then allows the proportion of adaptive substitutions (α) and the relative rate of adaptive substitution (ωa, expressed relative to the neutral substitution rate) to be estimated by combining the parameter estimates of the DFE with between-species nucleotide divergence data.
#	Under either the folded or unfolded SFS modes of operation, a simple model of recent demographic change is assumed, since demographic change can affect the selected and neutral SFSs in ways that resemble selection. A one-epoch unchanging population size, a two-epoch model with a change in population size from N1 to N2 t2 generations in the past, or a three-epoch model with two step changes from N1 to N2 and from N2 to N3 t2 and t3 generations in the past, respectively, can be run. Note that inferring the parameters of the 3-epoch model usually takes several hours of processing time.

#	DFE-alpha comprises three programs:#	1 - est_dfe – estimates the DFE for deleterious mutations. If analysing the unfolded SFS, the rate and fitness effects of advantageous mutation are also estimated.#	2 - est_alpha_omega – estimates the proportion of adaptive substitutions (α) and the relative rate of adaptive substitution (ωa) from the folded SFS.#	3 - prop_muts_in_s_ranges – estimates the proportions of deleterious mutations with fitness effects in different ranges of fitness effects on a scale Nes.

# the 3 different prompts with their command line options:
#	1 - ./est_dfe 				-c est_dfe_config_file.txt#	2 - ./est_alpha_omega 		-c est_alpha_omega_config_file.txt#	3 - ./prop_muts_in_s_ranges -c est_dfe_output_file 				-o output_file
#		i.e. this number 3 takes for input the output of program 1


#	The single input file for est_dfe contains one or more pairs of SFSs for selected and neutral sites, and has the following format (comments in red):# No. SFSs with different numbers of alleles sampled (m)	#	[integer on one line] for (i = 1 to m) {# No. alleles sampled in SFS i (xi) 						#	[= no. elements in unfolded vector]
# Selected SFS vector 										#	[a list of 0..xi counts]# Neutral SFS vector 										#	[a list of 0..xi counts]														#}
#	With the exception of the SFS vectors, all data elements are single integers on separate lines.
#	The SFS includes all sites and not just polymorphic sites. An SFS vector looks like this: number of sites with frequency 0/xi, # of sites with frequency 1/xi, ... number of sites with frequency xi/ xi
#	The SFS is a list of space-separated numbers (integer or real number) on a single line.
#	If folded SFSs are provided by the user, then their length must be the number of alleles sampled +1, and the upper half of the SFS vectors should be zeros.
#	Even if the neutral or selected SFS only are analysed in a given run, both SFSs need to be provided in the file



setwd("~/Documents/My_Documents/UofToronto/DFE_alpha/TestData")


generation <- 20000
polydat <- read.table("PolyMuts_Gen20000_outc_rep9.txt")
names(polydat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")

genome.dat <- read.table("Genomes_Gen20000_outc_rep9.txt", sep="A")

fixed.dat <- read.table("FixedData_outc_rep9.txt")
names(fixed.dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")


num.inds.sampled <- 100



	# because diploid:
	sample.size <- 2 * num.inds.sampled
	
	# make data frames of all possible neutral and deleterious mutations at all time points recorded
	neut.muts <- NULL
	seln.muts <- NULL
			
	## WHERE polydat IS A FILE OF ROWS OF MUTATIONS OCCURRING
	#	m1 = neutral site in coding
	#	m2, 3 = deleterious selected site in coding
	#	m4 = beneficial selected site in coding
		
# tack on fixed data and then can include counts for 0's
	fixed.mut.dat <- fixed.dat[fixed.dat$gen.fixed <= as.numeric(generation) ,]
		# this gives only mutations that have fixed PRIOR to and INCLUDING WITHIN the current generation time point sampled
	fixed.neut.muts <- c(which(fixed.mut.dat$mut.type == "m1"))
	fixed.seln.mut.IDs <- fixed.mut.dat$mut.ID[-fixed.neut.muts]
	fixed.neut.mut.IDs <- fixed.mut.dat$mut.ID[fixed.neut.muts]
	
	num.neut.muts.fixed <- length(fixed.neut.mut.IDs)
	num.seln.muts.fixed <- length(fixed.seln.mut.IDs)
	
	
	neut.muts <- polydat[polydat$mut.type == "m1" ,]
	seln.muts <- polydat[polydat$mut.type != "m1" ,]
			
	if(use.manual.sample == FALSE){
		sfs.total <- table(polydat$mut.prev)
		sfs.neut <- table(neut.muts$mut.prev)
		sfs.seln <- table(seln.muts$mut.prev)
		
		# add on fixed things to the fixed section of the table
		sfs.total[as.character(sample.size)] <- sfs.total[as.character(sample.size)] + num.neut.muts.fixed + num.seln.muts.fixed
		sfs.neut[as.character(sample.size)] <- sfs.neut[as.character(sample.size)] + num.neut.muts.fixed
		sfs.seln[as.character(sample.size)] <- sfs.seln[as.character(sample.size)] + num.seln.muts.fixed
		
		sfs.total["0"] <-
		sfs.total["0"] <-
		sfs.total["0"] <-
	
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





