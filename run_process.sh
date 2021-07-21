#!/usr/bin/env bash

# USER-SPECIFIED INPUTS, change these if necessary before running:

# Input directory where strain metadata and cluster assignment files are found
INPUTDIR=inputs

# Strain metadata file, with "Strain" identifier column
STRAINS=$INPUTDIR/strain_info.txt

# Unprocessed cluster assignments, with "Strain" identifier column
TP1_CLUSTERS=$INPUTDIR/tp1_clusters_init.txt
TP2_CLUSTERS=$INPUTDIR/tp2_clusters_init.txt
FORM=$INPUTDIR/form_inputs.txt

# Height of interest
HEIGHT=0

# Type of distances to save full set for
TYPEDIST="temp"

# Two pairs of coefficient sets, each in order "source, temporal geographic". 
# Example: "010" means source = 0, temporal = 1, geographic = 0
PARAMS="0-1-0,0-0-1"


# Do not change these variables unless necessary:
pdata=0


printf "\n\nPart 0/5:"
Rscript environment_setup.R


printf "\n\nPart 1/5:"
if [ -f $STRAINS -a -f $TP1_CLUSTERS -a -f $TP2_CLUSTERS -a -f $FORM ]; then
	Rscript scripts/1_prepare_inputs.R --inputdir $INPUTDIR --metadata $STRAINS --tp1 $TP1_CLUSTERS --tp2 $TP2_CLUSTERS --details $FORM
else
	echo "Not all required files found for part 1."
fi


if [ -f $INPUTDIR/processed/tp1_clusters.txt -a -f $INPUTDIR/processed/tp2_clusters.txt ]; then
	tp1_data=$INPUTDIR/processed/tp1_clusters.txt
	tp2_data=$INPUTDIR/processed/tp2_clusters.txt
	pdata=1
else
	echo "Data not processed properly"
fi


printf "\n\n\nPart 2/5:"
if [ $pdata == 1 ]; then
	Rscript scripts/2_cgm_collection.R --tp1 $tp1_data --tp2 $tp2_data --heights $HEIGHT
else
	echo "Not all required files found for part 2."
fi


printf "\n\n\nPart 3/5:"
if [ $pdata == 1 ]; then
	Rscript scripts/3_dist_matrices.R --strains $STRAINS --tp1 $tp1_data --tp2 $tp2_data --heights $HEIGHT --cpus 1 --trio $PARAMS
else
	echo "Not all required files found for part 3."
fi


printf "\n\n\nPart 4/5:"
if [ $pdata == 1 ]; then
	Rscript scripts/4_ecc_collection.R --strains $STRAINS --tp1 $tp1_data --tp2 $tp2_data --heights $HEIGHT --cpus 1 --trio $PARAMS
else
	echo "Not all required files found for part 4."
fi


printf "\n\n\nPart 5/5:"
if [ -f "results/ECCs.tsv" -a -f "results/CGM_strain_results.tsv" -a -f $STRAINS ]; then
  Rscript scripts/5_combine_results.R --ECCs results/ECCs.tsv --tp1 $INPUTDIR/processed/allTP1.Rds --tp2 $INPUTDIR/processed/allTP2.Rds --CGMs results/CGM_strain_results.tsv --strains $STRAINS --details $FORM
else
  echo "Not all required files found for part 5."
fi

# # Optional: note, this can only be run after 3_dist_matrices.R has been run successfully
# printf "\n\n\nOptional:"
# if [ -f "results/ECCs.tsv" -a -f "results/CGM_strain_results.tsv" -a -f $STRAINS ]; then
#   Rscript scripts/6_full_dists.R --strains $STRAINS --tp1 $tp1_data --tp2 $tp2_data --heights $HEIGHT --typedist $TYPEDIST
# else
#   echo "Not all required files found for part 6."
# fi

                            
printf "\nFinished process.\n"