
# make PBS scripts for specific westgrid servers given Slim input scripts


make.slim.input <- function(filename.start, rand.seed="1234567890", pop.size=10000, mut.rate, recomb.rate, mate.sys, prop.mate.type="", total.N.gens=10, samp.size=100, samp.type, rep){

	# options:
	#	random seed
	#	population size
	#	mutation rate
	#	recombination rate
	#	mating system
	#	optional - proportion mating type
	#	number of gens
	#	number inds to sample
	#	sample haploids or haploids
	
	
	sect1 <- "initialize() {
"
	sect2 <- paste(c("	setSeed(", rand.seed, ");
"), collapse="")
	sect3 <- paste(c("
	initializeMutationRate(", mut.rate, ");
"), collapse="")
	sect4 <- '
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", 0.3, "g", -0.01, 0.3);		// delet, mostly small effect
	initializeMutationType("m3", 0.05, "g", -0.5, 10);		// delet, few large effect

	initializeGenomicElementType("g1", c(m1,m2,m3), c(0.25, 0.7125, 0.0375 ));

	for (index in 0:99999){
		initializeGenomicElement(g1, index*1000, index*1000 + 199);
	}
'
	sect5 <- paste(c("
	initializeRecombinationRate(", recomb.rate, ");"), collapse="")
	sect6 <- '
}

1 {
'
	if(mate.sys == "outc"){
		sect7 <- paste(c('
	sim.addSubpop("p1", ', pop.size, ");"), collapse="")
	}
	if(mate.sys == "asex"){
		sect7 <- paste(c('
	sim.addSubpop("p1", ', pop.size, ");
	p1.setCloningRate(", prop.mate.type, ");"), collapse="")
	}
	if(mate.sys == "self"){
		sect7 <- paste(c('
	sim.addSubpop("p1", ', pop.size, ");
	p1.setSelfingRate(", prop.mate.type, ");"), collapse="")
	}
	sect8 <- "
}

"
	sampling.points <- format(seq(pop.size, (total.N.gens*pop.size), by= pop.size), scientific=FALSE)
	last.sample.point <- length(sampling.points)
	sect9 <- NULL
	for(i in 1:(last.sample.point-1)){
		if(samp.type == "haploid"){
			temp.sect <- paste(c(sampling.points[i], " late() { p1.outputSample(", samp.size, "); }"), collapse="")
		}
		if(samp.type == "diploid"){
			temp.sect <- paste(c(sampling.points[i], " late() { 
	subsampDiploids = sim.subpopulations.individuals;
	sampledIndividuals = sample(subsampDiploids, ", samp.size, ");
	sampledIndividuals.genomes.output();
}"), collapse="")
		}
		sect9 <- paste(c(sect9, temp.sect), collapse="\n")
	}
	
	sect10 <- "\n" 
	
	sect11 <- paste(c(sampling.points[last.sample.point], ' late() { sim.outputFull("FullOutput_', filename.start, "_N", pop.size, "_del_", mate.sys, prop.mate.type, "_rep", rep,'.txt"); }
',
	sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, "_N", pop.size, "_del_", mate.sys, prop.mate.type, "_rep", rep,'.txt"); }'), collapse="")
	
	

	file.text <- paste(c(sect1, sect2, sect3, sect4, sect5, sect6, sect7, sect8, sect9, sect10, sect11), collapse="")
	write(file.text, file=paste(c(filename.start, "_N", pop.size, "_del_", mate.sys, prop.mate.type, "_rep", rep,".txt"), collapse=""))
}	

