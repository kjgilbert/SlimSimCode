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
if echo $base_name | grep 20mbp; then
	genosize=20000000
if echo $base_name | grep 24mbp; then
	genosize=24000000
if echo $base_name | grep 26mbp; then
	genosize=26000000
Rscript CommandLine_RunSlimToDFEconversion.R $base_name subsample $dir genosize
Rscript CommandLine_RunSlimToDFEconversion.R $base_name full $dir genosize
# create divergence file (for all of them because easier to do in this loop anyway
Rscript CommandLine_RunSlimToAlphaOmega.R $base_name subsample $dir genosize
done

# now all the SFS inputs for DFE are ready
# go through each and do the DFE analyses

# put all those file names in a list
ls $dir | grep DFE_SFS > DFE_InputNames.txt


for input in `cat DFE_InputNames.txt`
do
# is the file we're on containing beneficials? if so, make the right type of config file
if echo $dir/$input | grep ben; then
	do_beneficial=TRUE
# make a config file
Rscript MakeConfigFiles_DFE.R $do_beneficial $dir $input -0.1 0.5
# run DFE
	# run class 0
./est_dfe -c $dir/config_class0.txt
	# run class 1
./est_dfe -c $dir/config_class1.txt
	# run Nes ranges
./prop_muts_in_s_ranges -c  results_dir_sel_class1/est_dfe.out -o output_file
# run DEF alpha_omega if doing that analysis
if test "$do_beneficial" = "true" ; then
	# run dfe_alpha_omega
	./est_alpha_omega -c $dir/config_alpha_omega.txt
done
