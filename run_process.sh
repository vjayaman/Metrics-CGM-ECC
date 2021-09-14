#!/usr/bin/env bash

# use sed -i 's/\r//g' run_process.sh to remove Windows carriage returns

# Strain metadata file, with "Strain" identifier column
RAWSTRAINS=inputs/strain_info.txt

# Unprocessed cluster assignments, with "Strain" identifier column
TP2_CLUSTERS=inputs/tp2_clusters_init.txt
FORM=inputs/form_inputs.txt

# DO NOT CHANGE pdata, it's a flag used to check for success
pdata=0

printf "\n\nPart 0/6:"
Rscript environment_setup.R


printf "\n\nPart 1/6:"
if [ -f $RAWSTRAINS -a -f $TP2_CLUSTERS -a -f $FORM ]; then
	Rscript scripts/1_prepare_inputs.R -m $RAWSTRAINS -b $TP2_CLUSTERS -d $FORM
else
	echo "Not all required files found for part 1."
fi


if [ -f inputs/processed/tp2_clusters.txt -a -f inputs/processed/clustersets.Rds ]; then
	tp2_data=inputs/processed/tp2_clusters.txt
	interval_file=inputs/processed/clustersets.Rds
	STRAINS=inputs/processed/strain_info.txt
	pdata=1
else
	echo "Data not processed properly"
fi


printf "\n\n\nPart 2/6:"
if [ $pdata == 1 ]; then
	Rscript scripts/2_cgm_collection.R -f $interval_file -d $FORM
else
	echo "Not all required files found for part 2."
fi


printf "\n\n\nPart 3/6:"
if [ $pdata == 1 ]; then
	Rscript scripts/3_dist_matrices.R -m $STRAINS -b $tp2_data -f $interval_file -d $FORM
else
	echo "Not all required files found for part 3."
fi


printf "\n\n\nPart 4/6:"
if [ $pdata == 1 ]; then
	Rscript scripts/4_averages.R -m $STRAINS -b $tp2_data -f $interval_file -d $FORM
else
	echo "Not all required files found for part 4."
fi


printf "\n\n\nPart 5/6:"
if [ $pdata == 1 ]; then
	Rscript scripts/5_ecc_collection.R -m $STRAINS -b $tp2_data -f $interval_file -d $FORM
else
	echo "Not all required files found for part 4."
fi


printf "\n\n\nPart 6/6:"
# if [ -f "results/ECCs.tsv" -a -f "results/CGM_strain_results.tsv" -a -f $STRAINS ]; then
  Rscript scripts/6_combine_results.R -m $STRAINS -f $interval_file -d $FORM
# else
#   echo "Not all required files found for part 5."
# fi
              
# printf "\n\n\nPart 6/6:"
# if [ -f "results/CGM_strain_results.tsv" -a -f "intermediate_data/dist_extremes.Rds" -a -f $STRAINS ]; then
# 	Rscript scripts/6_heatmap_dists.R -m $STRAINS -b $tp2_data -c results/CGM_strain_results.tsv -d $FORM
# else
# 	echo "Not all required files found for part 6."
# fi

                            
printf "\nFinished process.\n"