source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/MakePBSfiles_WestgridSlim.R', chdir = TRUE)
#source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/MakeInputFiles_Slim.R', chdir = TRUE)
# to instead make it with beneficial mutations as well:
source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/MakeInputFiles_Slim_specifyGenomeBens_TransMateSysFROMexistingOutput.R', chdir = TRUE)



dir <- "~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Nov3_TransitionMatingSystems/Inputs_Nov2_N10000"

setwd(dir)

existing.dir <- "~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Nov3_TransitionMatingSystems/Inputs_Nov2_N10000/ExistingOutputs"
	# this is the directory on my computer with the full output files that it will use the name of to write the new inputs - it must only contain the outcrossing full output files in the dir

existing.outputs <- system(paste(c("ls ", existing.dir), collapse=""), intern=TRUE)


rand.seeds <- c(
"1234567890",
"2345678901",
"3456789012",
"4567890123",
"5678901234",
"6789012345",
"7890123456",
"8901234567",
"9012345678",
"0123456789"
)

reps <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

mut.rate <- "7e-9"
recomb.rate <- "5e-8"
sample.dips <- "diploid"
N.number.generations <- 10
inds.to.sample <- 100


for(i in existing.outputs){
	filename <- unlist(strsplit(i, split="_"))
	
	replicate <- unlist(strsplit(filename[grep("rep", filename)], split="\\."))[1]
	if(replicate == "rep1"){ rand.sd <- rand.seeds[1]; repl <- reps[1]}
	if(replicate == "rep2"){ rand.sd <- rand.seeds[2]; repl <- reps[2] }
	if(replicate == "rep3"){ rand.sd <- rand.seeds[3]; repl <- reps[3] }
	if(replicate == "rep4"){ rand.sd <- rand.seeds[4]; repl <- reps[4] }
	if(replicate == "rep5"){ rand.sd <- rand.seeds[5]; repl <- reps[5] }
	if(replicate == "rep6"){ rand.sd <- rand.seeds[6]; repl <- reps[6] }
	if(replicate == "rep7"){ rand.sd <- rand.seeds[7]; repl <- reps[7] }
	if(replicate == "rep8"){ rand.sd <- rand.seeds[8]; repl <- reps[8] }
	if(replicate == "rep9"){ rand.sd <- rand.seeds[9]; repl <- reps[9] }
	if(replicate == "rep10"){ rand.sd <- rand.seeds[10]; repl <- reps[10] }
	
	if(filename[3] == "N10000") pop.size <- 10000
	if(filename[4] == "30mbp") genosize <- "30mbp"

	if(filename[5] == "ben-del"){ bens <- TRUE }else{ bens <- FALSE}

	westgrid.stored.dir <- "/home/kgilbert/Slim/Outputs/Done_Oct8_Outputs"
		# this is the directory on westgrid where it will find the file when running slim from that starting point
	
	make.slim.input(filename.start="Nov3", 
		rand.seed= rand.sd, pop.size=pop.size, genome.size=genosize, mut.rate=mut.rate, ben.muts=bens, recomb.rate=recomb.rate, 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=repl,
		mate.sys="self", prop.mate.type="0.90", input.file=paste(c(westgrid.stored.dir, "/", i), collapse="")
	)
		make.slim.input(filename.start="Nov3", 
		rand.seed= rand.sd, pop.size=pop.size, genome.size=genosize, mut.rate=mut.rate, ben.muts=bens, recomb.rate=recomb.rate, 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=repl,
		mate.sys="self", prop.mate.type="0.99", input.file=paste(c(westgrid.stored.dir, "/", i), collapse="")
	)
	westgrid.stored.dir <- "/global/scratch/kgilbert/Slim/Inputs/Temp"
		# this is the directory on westgrid where it will find the file when running slim from that starting point

	make.slim.input(filename.start="Nov3", 
		rand.seed= rand.sd, pop.size=pop.size, genome.size=genosize, mut.rate=mut.rate, ben.muts=bens, recomb.rate=recomb.rate, 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=repl,
		mate.sys="asex", prop.mate.type="0.90", input.file=paste(c(westgrid.stored.dir, "/", i), collapse="")
	)
		make.slim.input(filename.start="Nov3", 
		rand.seed= rand.sd, pop.size=pop.size, genome.size=genosize, mut.rate=mut.rate, ben.muts=bens, recomb.rate=recomb.rate, 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=repl,
		mate.sys="asex", prop.mate.type="0.99", input.file=paste(c(westgrid.stored.dir, "/", i), collapse="")
	)
}





# see if these will work on grex:

##	RAM <- "18gb"
##	walltime <- "168:00:00"
##	server <- "grex"

##	multi.pbs(input.file.directory=dir, RAM=RAM, walltime=walltime, server=server)



# do the faster ones on jasper:

ini.dir <- "~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Nov3_TransitionMatingSystems/Inputs_Nov2_N10000"

RAM <- "6gb"
walltime <- "72:00:00"
server <- "jasper"
multi.pbs(input.file.directory=ini.dir, RAM=RAM, walltime=walltime, server=server)


ini.dir <- "~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Nov3_TransitionMatingSystems/Inputs_Nov2_N10000/asex99_grex"

RAM <- "14gb"
walltime <- "165:00:00"
server <- "grex"
multi.pbs(input.file.directory=ini.dir, RAM=RAM, walltime=walltime, server=server)





