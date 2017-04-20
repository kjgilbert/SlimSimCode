
make.slim.input <- function(filename.start, input.file, rand.seed="1234567890", pop.size=10000, ben.muts=FALSE, delet.effect, ben.effect, mut.props, mate.sys, prop.mate.type="", scenario.name, total.N.gens=10, samp.size=100, samp.type, sampling.points, rep, dontSampleFull=FALSE, dominance=0.3){
	
	sect1 <- 'initialize() {
'
	sect2 <- paste(c(' setSeed(', rand.seed, ');
'), collapse="")
	sect3 <- paste(c('
	initializeMutationRate(2.5e-8);
'), collapse="")


if(ben.muts == FALSE){
	sect4 <- paste(c('
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", ', dominance, ', "g", ', delet.effect, ', 0.2);		// delet, mostly small effect

// if want 75% of mutations to be deleterious:
	initializeGenomicElementType("g1", c(m1,m2), c(0.25, 0.75));

'), collapse="")
}
if(ben.muts == TRUE){
	sect4 <- paste(c('
	initializeMutationType("m1", 0.5, "f", 0.0);			// neutral 
	initializeMutationType("m2", ', dominance, ', "g", ', delet.effect, ', 0.2);		// delet, mostly small effect
	initializeMutationType("m3", 0.5, "g", ', ben.effect, ', 0.2);		// beneficial

// if want 75% of mutations to be deleterious:
	initializeGenomicElementType("g1", c(m1,m2,m3), c(0.25, ', mut.props, '));

'), collapse="")
}
     
sect5 <- paste("
        
        for ( index in 0:249){
		initializeGenomicElement(g1, index*40000, index*40000+700);
		initializeGenomicElement(g1, index*40000+2200, index*40000+2350);
		initializeGenomicElement(g1, index*40000+3850, index*40000+4000);
		initializeGenomicElement(g1, index*40000+5500, index*40000+5650);
		initializeGenomicElement(g1, index*40000+7150, index*40000+7300);
		initializeGenomicElement(g1, index*40000+8800, index*40000+8950);
		initializeGenomicElement(g1, index*40000+10450, index*40000+10600);
		initializeGenomicElement(g1, index*40000+12100, index*40000+12500);
	}
    initializeRecombinationRate(1e-8, 9999999);
}

1 {
")

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
	## input the sampling points instead		#sampling.points <- format(seq(pop.size, (total.N.gens*pop.size), by= pop.size), scientific=FALSE)
	last.sample.point <- length(sampling.points)
	sect9 <- NULL
	for(i in 1:(last.sample.point-1)){
		if(samp.type == "haploid"){
			temp.sect <- paste(c(sampling.points[i], " late() { p1.outputSample(", samp.size, "); }"), collapse="")
		}
		if(samp.type == "diploid"){
			temp.sect <- paste(c(sampling.points[i], ' late() { 
	subsampDiploids = sim.subpopulations.individuals;
	sampledIndividuals = sample(subsampDiploids, ', samp.size, ');
	sampledIndividuals.genomes.output();
}'), collapse="")
		}
		sect9 <- paste(c(sect9, temp.sect), collapse="\n")
	}
	
	sect10 <- '\n'
	
if(ben.muts == FALSE){
	sect11 <- paste(c(sampling.points[last.sample.point], ' late() { sim.outputFull("FullOutput_', filename.start, '_N', pop.size, '_10mbp_del_', scenario.name, mate.sys, prop.mate.type, '_rep', rep, '.txt"); }
',
	sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_10mbp_del_', scenario.name, mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")

if(dontSampleFull==TRUE & samp.type == 'diploid'){
sect11 <- paste(c(sampling.points[last.sample.point], ' late() { 
	subsampDiploids = sim.subpopulations.individuals;
	sampledIndividuals = sample(subsampDiploids, ', samp.size*10, ');
	sampledIndividuals.genomes.output();
}
',
	sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_10mbp_del_', scenario.name, mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")
}

	file.text <- paste(c(sect1, sect2, sect3, sect4, sect5, sect6, sect7, sect7p5, sect8, sect9, sect10, sect11), collapse="")
	write(file.text, file=paste(c(filename.start, "_N", pop.size, "_10mbp_del_", scenario.name, mate.sys, prop.mate.type, "_rep", rep,".txt"), collapse=""))
}
	
if(ben.muts == TRUE){
	sect11 <- paste(c(sampling.points[last.sample.point], ' late() { sim.outputFull("FullOutput_', filename.start, '_N', pop.size, '_10mbp_ben-del_', scenario.name, mate.sys, prop.mate.type, '_rep', rep, '.txt"); }
',
	sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_10mbp_ben-del_', scenario.name, mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")

if(dontSampleFull==TRUE & samp.type == "diploid"){
sect11 <- paste(c(sampling.points[last.sample.point], ' late() {
 	subsampDiploids = sim.subpopulations.individuals;
	sampledIndividuals = sample(subsampDiploids, ', samp.size*10, ');
	sampledIndividuals.genomes.output();
}
',
	sampling.points[last.sample.point], ' late() { sim.outputFixedMutations("FixedOutput_', filename.start, '_N', pop.size, '_10mbp_ben-del_', scenario.name, mate.sys, prop.mate.type, '_rep', rep, '.txt"); }'), collapse="")
}

	file.text <- paste(c(sect1, sect2, sect3, sect4, sect5, sect6, sect7, sect7p5, sect8, sect9, sect10, sect11), collapse="")
	write(file.text, file=paste(c(filename.start, "_N", pop.size, "_10mbp_ben-del_", scenario.name, mate.sys, prop.mate.type, "_rep", rep,".txt"), collapse=""))
}		
	
}



#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>




dir <- "~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/April2017_NewMesserLikeSims/Transitions_10N-25N/plasmodiumBashScripts_AndInputs"

setwd(dir)

rand.seeds <- c("1234567890","2345678901","3456789012","4567890123","5678901234","6789012345","7890123456", "8901234567", "9012345678", "0123456789")
reps <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")


existing.outputs <- c(
	"FullOutput_Apr19_10N_b01d01_N10000_10mbp_ben-del_outc_rep", 
	"FullOutput_Apr19_10N_b01d0005_N10000_10mbp_ben-del_outc_rep", 
	"FullOutput_Apr19_10N_b001d01_N10000_10mbp_ben-del_outc_rep", 
	"FullOutput_Apr19_10N_b001d0005_N10000_10mbp_ben-del_outc_rep", 
	"FullOutput_Apr19_10N_d01_N10000_10mbp_del_outc_rep", 
	"FullOutput_Apr19_10N_d0005_N10000_10mbp_del_outc_rep")

# their.directory <- "/fast_disk/kgilbert/Inputs/NewMesser_TransitionStartFiles/"
their.directory <- "/plas1/kgilbert/Slim/Inputs/NewMesser_TransitionStartFiles/"

pop.size <- 10000
sample.dips <- "diploid"
N.number.generations <- 15
inds.to.sample <- 100

## partial dominance
dom.val <- 0.3

samp.pts <- format(c((pop.size*10)+1, seq((pop.size*10.01),(pop.size*10.2), by=100), seq((pop.size*10.3), (pop.size*10.5), by=1000), seq((pop.size*11), (pop.size*25), by=10000)), scientific=FALSE)



for(i in existing.outputs){
	iterate <- 1
	
	for(j in reps){
		filename <- paste(c(i, j, ".txt"), collapse="")
		rand.sd <- rand.seeds[iterate]
		
		if(length(grep("_b001d0005", filename)) == 1) { mut.names <- "b001d0005" ; s.d <- "-0.0005" ; s.b="0.001"; mut.prop <- "0.749375, 0.000625" }
		if(length(grep("_b001d01", filename)) == 1) { mut.names <- "b001d01" ; s.d <- "-0.01" ; s.b="0.001"; mut.prop <- "0.749375, 0.000625" }
		if(length(grep("_b01d0005", filename)) == 1) { mut.names <- "b01d0005" ; s.d <- "-0.0005" ; s.b="0.01"; mut.prop <- "0.749937, 0.000063" }
		if(length(grep("_b01d01", filename)) == 1) { mut.names <- "b01d01" ; s.d <- "-0.01" ; s.b="0.01"; mut.prop <- "0.749937, 0.000063" }
		if(length(grep("_d0005", filename)) == 1) { mut.names <- "d0005" ; s.d <- "-0.0005" ; s.b=NULL }
		if(length(grep("_d01", filename)) == 1) { mut.names <- "d01" ; s.d <- "-0.01" ; s.b=NULL }
	
		if(is.null(s.b)){
			#outcrossers
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=FALSE, mate.sys="outc", 
				scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				dominance=dom.val, delet.effect=s.d, sampling.points=samp.pts)	
		
			# selfers
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=FALSE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="self", prop.mate.type="0.90", dominance=dom.val, delet.effect=s.d, sampling.points=samp.pts)	
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=FALSE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="self", prop.mate.type="0.99", dominance=dom.val, delet.effect=s.d, sampling.points=samp.pts)	
			#asex
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=FALSE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="asex", prop.mate.type="0.90", dominance=dom.val, delet.effect=s.d, sampling.points=samp.pts)	
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=FALSE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="asex", prop.mate.type="0.99", dominance=dom.val, delet.effect=s.d, sampling.points=samp.pts)	
		}else{
			#outcrossers
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=TRUE, mate.sys="outc", 
				scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				dominance=dom.val, delet.effect=s.d, ben.effect=s.b, mut.props=mut.prop, sampling.points=samp.pts)
			#selfers
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=TRUE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="self", prop.mate.type="0.90", dominance=dom.val, delet.effect=s.d, ben.effect=s.b, mut.props=mut.prop, sampling.points=samp.pts)
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=TRUE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="self", prop.mate.type="0.99", dominance=dom.val, delet.effect=s.d, ben.effect=s.b, mut.props=mut.prop, sampling.points=samp.pts)
			#asex
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=TRUE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="asex", prop.mate.type="0.90", dominance=dom.val, delet.effect=s.d, ben.effect=s.b, mut.props=mut.prop, sampling.points=samp.pts)
			make.slim.input(filename.start=paste(c("Apr21_25N_", mut.names), collapse=""), input.file=paste(c(their.directory, filename), collapse=""),
				rand.seed=rand.sd, pop.size=pop.size, ben.muts=TRUE, scenario.name="Transition-",
				total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate], 
				mate.sys="asex", prop.mate.type="0.99", dominance=dom.val, delet.effect=s.d, ben.effect=s.b, mut.props=mut.prop, sampling.points=samp.pts)

		}
		iterate <- iterate + 1
	}
}


