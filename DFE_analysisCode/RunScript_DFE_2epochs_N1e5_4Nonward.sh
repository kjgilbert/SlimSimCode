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

# put all those file names in a list
ls $dir | grep Fixed | sed "s/FixedOutput_//g" > base_InputNames.txt

# loop through the list and do each one at a time
for base_name in `cat base_InputNames.txt`
do
if echo $base_name | grep 20mbp
then
	genosize=20000000
fi
if echo $base_name | grep 24mbp
then
	genosize=24000000
fi
if echo $base_name | grep 25mbp
then
	genosize=25000000
fi
if echo $base_name | grep 26mbp
then
	genosize=26000000
fi
if echo $base_name | grep 30mbp
then
	genosize=30000000
fi
# turn the sample files into a format suitable for a "FullOutput..." file
Rscript CommandLine_MakeModifiedFullSampleFile_ForDFEinput.R $base_name $dir
# then analyze like normal:
Rscript CommandLine_RunSlimToDFEconversion_N1e5_4Nonward.R $base_name subsample $dir $genosize
# create divergence file (for all of them because easier to do in this loop anyway
Rscript CommandLine_RunSlimToAlphaOmega_N1e5_4Nonward.R $base_name subsample $dir $genosize
done

# now all the SFS inputs for DFE are ready
# go through each and do the DFE analyses

# put all those file names in a list
ls $dir | grep DFE_SFS_sub > DFE_InputNames.txt


for input in `cat DFE_InputNames.txt`
do
# is the file we're on containing beneficials? if so, make the right type of config file
##	just do beneficial estimation the whole time
##	if echo $dir/$input | grep ben
##	then
	do_beneficial=TRUE
##	fi
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
