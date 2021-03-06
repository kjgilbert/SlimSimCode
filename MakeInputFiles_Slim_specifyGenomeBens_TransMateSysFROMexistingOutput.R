
# make PBS scripts for specific westgrid servers given Slim input scripts


make.slim.input <- function(filename.start, rand.seed="1234567890", pop.size=10000, genome.size, mut.rate, ben.muts=FALSE, recomb.rate, mate.sys, prop.mate.type="", total.N.gens=10, sampling.points, samp.size=100, samp.type, rep, dontSampleFull=FALSE, newS=FALSE, input.file){

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


if(ben.muts == FALSE){
	sect4 <- '
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", 0.3, "g", -0.01, 0.3);		// delet, mostly small effect
	initializeMutationType("m3", 0.05, "g", -0.5, 10);		// delet, few large effect

// if want 75% of mutations to be deleterious:
	initializeGenomicElementType("g1", c(m1,m2,m3), c(0.25, 0.7125, 0.0375 ));

'
	if(newS == TRUE){
			sect4 <- '
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", 0.3, "g", -0.1, 0.3);		// delet, mostly small effect
	initializeMutationType("m3", 0.05, "g", -5, 10);		// delet, few large effect

// if want 75% of mutations to be deleterious:
	initializeGenomicElementType("g1", c(m1,m2,m3), c(0.25, 0.7125, 0.0375 ));

'
	}
}
if(ben.muts == TRUE){
	sect4 <- '
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", 0.3, "g", -0.01, 0.3);		// delet, mostly small effect
	initializeMutationType("m3", 0.05, "g", -0.5, 10);		// delet, few large effect
	initializeMutationType("m4", 0.5, "g", 0.01, 0.3);		// beneficial

// if want 75% of mutations to be deleterious:
	initializeGenomicElementType("g1", c(m1,m2,m3,m4), c(0.25, 0.7125, 0.0375, 0.00071 ));

'
	if(newS == TRUE){
	sect4 <- '
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", 0.3, "g", -0.1, 0.3);		// delet, mostly small effect
	initializeMutationType("m3", 0.05, "g", -5, 10);		// delet, few large effect
	initializeMutationType("m4", 0.5, "g", 0.1, 0.3);		// beneficial

// if want 75% of mutations to be deleterious:
	initializeGenomicElementType("g1", c(m1,m2,m3,m4), c(0.25, 0.7125, 0.0375, 0.00071 ));

'
	}
}
if(genome.size == "20mbp"){
	sect5 <- paste(c("
        // one chromosome with coding elements over 200 bp then replace 800 bp noncoding by 800x recombination rate, up to a total of 20Mbp in size
        initializeGenomicElement(g1, 0, 19999999);

        for (index in 0:99999){
                initializeRecombinationRate(", recomb.rate, ", index*200);
                initializeRecombinationRate((800*(", recomb.rate, ")), index*200 + 1);
        }
        initializeRecombinationRate(", recomb.rate, ", 19999999);"), collapse="")
}  
if(genome.size == "24mbp"){
	sect5 <- paste(c("
        // one chromosome with coding elements over 200 bp then replace 800 bp noncoding by 800x recombination rate, up to a total of 20Mbp in size
        initializeGenomicElement(g1, 0, 23999999);

        for (index in 0:119999){
                initializeRecombinationRate(", recomb.rate, ", index*200);
                initializeRecombinationRate((800*(", recomb.rate, ")), index*200 + 1);
        }
        initializeRecombinationRate(", recomb.rate, ", 23999999);"), collapse="")
}      
if(genome.size == "25mbp"){
	sect5 <- paste(c("
        // one chromosome with coding elements over 200 bp then replace 800 bp noncoding by 800x recombination rate, up to a total of 20Mbp in size
        initializeGenomicElement(g1, 0, 24999999);

        for (index in 0:124999){
                initializeRecombinationRate(", recomb.rate, ", index*200);
                initializeRecombinationRate((800*(", recomb.rate, ")), index*200 + 1);
        }
        initializeRecombinationRate(", recomb.rate, ", 24999999);"), collapse="")
}
if(genome.size == "26mbp"){
	sect5 <- paste(c("
        // one chromosome with coding elements over 200 bp then replace 800 bp noncoding by 800x recombination rate, up to a total of 20Mbp in size
        initializeGenomicElement(g1, 0, 25999999);

        for (index in 0:129999){
                initializeRecombinationRate(", recomb.rate, ", index*200);
                initializeRecombinationRate((800*(", recomb.rate, ")), index*200 + 1);
        }
        initializeRecombinationRate(", recomb.rate, ", 25999999);"), collapse="")
}
if(genome.size == "30mbp"){
	sect5 <- paste(c('
        // one chromosome with coding elements over 200 bp then replace 800 bp noncoding by 800x recombination rate, up to a total of 20Mbp in size
        initializeGenomicElement(g1, 0, 29999999);

        for (index in 0:149999){
                initializeRecombinationRate(', recomb.rate, ', index*200);
                initializeRecombinationRate((800*(', recomb.rate, ')), index*200 + 1);
        }
        initializeRecombinationRate(', recomb.rate, ', 29999999);'), collapse="")
}      
	sect6 <- paste(c('
}

', format(pop.size*10, scientific=FALSE), ' late() {
	sim.readFromPopulationFile("', input.file, '");'), collapse="")
	
	sect7 <- paste(c('
}

', format((pop.size*10)+1, scientific=FALSE), ' {
'), collapse="")	

	if(mate.sys == "outc"){
		sect7p5 <- ""
	}
	if(mate.sys == "asex"){
		sect7p5 <- paste(c('
	p1.setCloningRate(', prop.mate.type, ");"), collapse="")
	}
	if(mate.sys == "self"){
		sect7p5 <- paste(c('
	p1.setSelfingRate(', prop.mate.type, ");"), collapse="")
	}
	sect8 <- "
}

"
##	sampling.points <- format(seq(((total.N.gens+1)*pop.size), (2*total.N.gens*pop.size), by= pop.size), scientific=FALSE)
	last.sample.point <- length(sampling.points)
	sect9 <- NULL
	for(i in 1:last.sample.point){
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
	
if(ben.muts == FALSE){
	sect11 <- paste(c(sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_', genome.size, '_del_TransTo', mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")

if(dontSampleFull==TRUE & samp.type == "diploid"){
sect11 <- paste(c(sampling.points[last.sample.point], " late() { 
	subsampDiploids = sim.subpopulations.individuals;
	sampledIndividuals = sample(subsampDiploids, ", samp.size*10, ");
	sampledIndividuals.genomes.output();
}
",
	sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_', genome.size, '_del_TransTo', mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")
}

	file.text <- paste(c(sect1, sect2, sect3, sect4, sect5, sect6, sect7, sect7p5, sect8, sect9, sect10, sect11), collapse="")
	write(file.text, file=paste(c(filename.start, "_N", pop.size, "_", genome.size, "_del_TransTo", mate.sys, prop.mate.type, "_rep", rep,".txt"), collapse=""))
}
	
if(ben.muts == TRUE){
	sect11 <- paste(c(sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_', genome.size, '_ben-del_TransTo', mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")

if(dontSampleFull==TRUE & samp.type == "diploid"){
sect11 <- paste(c(sampling.points[last.sample.point], " late() {
 	subsampDiploids = sim.subpopulations.individuals;
	sampledIndividuals = sample(subsampDiploids, ", samp.size*10, ");
	sampledIndividuals.genomes.output();
}
",
	sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_', genome.size, '_ben-del_TransTo', mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")
}

	file.text <- paste(c(sect1, sect2, sect3, sect4, sect5, sect6, sect7, sect7p5, sect8, sect9, sect10, sect11), collapse="")
	write(file.text, file=paste(c(filename.start, "_N", pop.size, "_", genome.size, "_ben-del_TransTo", mate.sys, prop.mate.type, "_rep", rep,".txt"), collapse=""))
}		
	

}	

