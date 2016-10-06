#!/bin/bash

# run DFE-alpha on Slim outputs

# for 1 individual slim output, need to:
#	calculate neutral and selected SFS
#	make config file
#	run est_dfe class 0
#	run est_dfe class 1 with output from class 0
#	run prop_muts_in_s_range on class 1 output

# if beneficials, run alpha_omega on est_dfe output class 0 and 1


# get all outputs in dir to run analyses on
echo "Provide the path to the directory containing the SLiM outputs to be analysed."
read dir



# now all the SFS inputs for DFE are ready
# go through each and do the DFE analyses

# put all those file names in a list
ls $dir | grep DFE_SFS_sub > DFE_InputNames.txt


for input in `cat DFE_InputNames.txt`
do
# is the file we're on containing beneficials? if so, make the right type of config file
##if echo $dir/$input | grep ben
##then
	do_beneficial=TRUE
##fi
# make the directories for the outputs from this analysis
basedir=$( echo $input | sed "s/.txt//g" )
basedir0="Outputs/OutputClass0_"
basedir1="Outputs/OutputClass1_"
dir0=${basedir0}${basedir}
dir1=${basedir1}${basedir}
mkdir $dir0
mkdir $dir1
# make a config file
Rscript MakeConfigFiles_DFE_2epochs.R $do_beneficial $dir $input -0.1 0.5
# run DFE
	# run class 0
./est_dfe -c ${dir}config_class0.txt
	# run class 1
./est_dfe -c ${dir}config_class1.txt
	# run Nes ranges
results_file_sel_class1=${basedir1}${basedir}"/est_dfe.out"
output_file="output_prop_muts_"${input}
./prop_muts_in_s_ranges -c $results_file_sel_class1 -o $output_file
# run DEF alpha_omega if doing that analysis
##	if test "$do_beneficial" = "true"
##	then
	# run dfe_alpha_omega
./est_alpha_omega -c $dir/config_alpha_omega.txt
##	fi
done
