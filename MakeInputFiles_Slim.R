
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
	
	
	section1 <- "initialize() {"
	section2 <- paste(c("	setSeed(", rand.seed, ");
		"), collapse="")
	section3 <- paste(c("	initializeMutationRate(", mut.rate, ");"), collapse="")
	section4 <- '
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", 0.3, "g", -0.01, 0.3);		// delet, mostly small effect
	initializeMutationType("m3", 0.05, "g", -0.5, 10);		// delet, few large effect

	initializeGenomicElementType("g1", c(m1,m2,m3), c(0.25, 0.7125, 0.0375 ));

	for (index in 0:99999){
		initializeGenomicElement(g1, index*1000, index*1000 + 199);
	}
'
	section5 <- paste(c("initializeRecombinationRate(", recomb.rate, ");"), collapse="")
	section6 <- '
}

1 {
'
	if(mate.sys == "outc"){
		section7 <- paste(c('sim.addSubpop("p1", ', pop.size, ");"), collapse="")
	}
	if(mate.sys == "asex"){
		section7 <- paste(c('sim.addSubpop("p1", ', pop.size, ");
	p1.setCloningRate(", prop.mate.type, ");"), collapse="")
	}
	if(mate.sys == "outc"){
		section7 <- paste(c('sim.addSubpop("p1", ', pop.size, ");
	p1.setSelfingRate(", prop.mate.type, ");"), collapse="")
	}
	section8 <- "
}

"
	
	
	sampling.points <- format(seq(pop.size, (total.N.gens*pop.size), by= pop.size), scientific=FALSE)


	file.text <- paste(c(section1, section2, section3, section4, section5, section6, section7), collapse="")
	write(file.text, file=paste(c(filename.start, "_N", pop.size, "_del_", mate.sys, prop.mate.type, "_rep", rep,".txt"), collapse=""))
}	


