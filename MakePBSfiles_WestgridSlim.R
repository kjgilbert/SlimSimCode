
# make PBS scripts for specific westgrid servers given Slim input scripts


make.pbs <- function(filename, RAM, walltime=NULL, westgrid.server=NULL){
	# currently have slim 2.1 installed on: orcinus, bugaboo, grex and jasper
	# orcinus, walltime limit = 240:00:00  (10 days)
	# bugaboo, walltime limit = 122 days, so as much as needed
	# grex, walltime limit = 168:00:00  (7 days)
	# jasper, walltime limit = 72:00:00  (3 days)

	if(is.null(westgrid.server)){
		print("Must specify 'grex', 'bugaboo', 'jasper', or 'orcinus' as the westgrid server being used.")
	}

	if(westgrid.server=="grex"){
		slim.directory <- "/global/scratch/kgilbert/Slim/Slim_v2.1/slim"
		input.directory <- "/global/scratch/kgilbert/Slim/Inputs/"
		output.directory <- "> /global/scratch/kgilbert/Slim/Outputs/"
		load.command <- "module load gcc/5.2.0-experimental"
	}
	if(westgrid.server=="orcinus"){
		slim.directory <- "/global/scratch/kgilbert/Slim/Slim_v2.1/bin/slim"
		input.directory <- "/global/scratch/kgilbert/Slim/Inputs/"
		output.directory <- "> /global/scratch/kgilbert/Slim/Outputs/"
		load.command <- "module load gcc/4.8.2rev203690"
	}
	if(westgrid.server=="bugaboo"){
		slim.directory <- "/global/scratch/kgilbert/Slim/Slim_v2.1/bin/slim"
		input.directory <- "/global/scratch/kgilbert/Slim/Inputs/"
		output.directory <- "> /global/scratch/kgilbert/Slim/Outputs/"
		load.command <- "module load gcc/4.8.0"
	}
	if(westgrid.server=="jasper"){
		slim.directory <- "/home/kgilbert/Slim/slim_v2.1"
		input.directory <- "/home/kgilbert/Slim/Inputs/"
		output.directory <- "> /home/kgilbert/Slim/Outputs/"
		load.command <- "module load compiler/gcc/4.8.2"
	}
		
	section1 <- "#!/bin/bash
#PBS -S /bin/bash
#PBS -l procs=1
#PBS -m bea
#PBS -M kgilbert@zoology.ubc.ca
"
	section2 <- paste(c("#PBS -l walltime=", walltime), collapse="")	
	section3 <- paste(c("#PBS -l pmem=", RAM), collapse="")	
	section4 <- '

cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"

echo "Starting run at: `date`"

'
	section5 <- paste(load.command)
	section5.5 <- "
"
	section6 <- paste(c(slim.directory, paste(c(input.directory, filename), collapse=""), paste(c(output.directory, "SampleOutput_", filename), collapse="")), collapse=" ")
	section7 <- '
echo "Job finished with exit code $? at: `date`"

'

	file.text <- paste(c(section1, section2, section5.5, section3, section4, section5, section5.5, section6, section7), collapse="")
	write(file.text, file=paste(c("PBS_", westgrid.server, "_", filename, ".pbs"), collapse=""))
}	


multi.pbs <- function(input.file.directory, RAM, walltime=NULL, server){
	
	setwd(input.file.directory)
	files.in.folder <- system("ls", intern=TRUE)
	# take all files in the directory
	
	for(i in 1:length(files.in.folder)){
		make.pbs(filename=files.in.folder[i], 
		RAM=RAM, 
		walltime=walltime, 
		westgrid.server=server)
	}
}

