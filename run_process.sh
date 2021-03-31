#!/usr/bin/env bash

# USER-SPECIFIED INPUTS, change these if necessary before running:

# Input directory where strain metadata and cluster assignment files are found
INPUTDIR=inputs

# Strain metadata file, with "Strain" identifier column
STRAINS=$INPUTDIR/strain_info.txt

# Unprocessed cluster assignments, with "Strain" identifier column
TP1_CLUSTERS=$INPUTDIR/tp1_clusters_init.txt
TP2_CLUSTERS=$INPUTDIR/tp2_clusters_init.txt

# Height of interest
HEIGHT=0

# Two pairs of coefficient sets, each in order "source, temporal geographic". 
# Example: "010" means source = 0, temporal = 1, geographic = 0
PARAMS="010-001" 


# Do not change these variables unless necessary:
pdata=0
SOURCE=use_placeholder


printf "\n\nPart 1/4:"
if [ -f $STRAINS -a -f $TP1_CLUSTERS -a -f $TP2_CLUSTERS ]; then
	Rscript scripts/prepare_inputs.R --inputdir $INPUTDIR --metadata $STRAINS --tp1 $TP1_CLUSTERS --tp2 $TP2_CLUSTERS --source $SOURCE
else
	echo "Not all required files found for part 1."
fi


if [ -f $INPUTDIR/processed/tp1_clusters.txt -a -f $INPUTDIR/processed/tp2_clusters.txt -a -f $INPUTDIR/processed/source_data.tsv ]; then
	tp1_data=$INPUTDIR/processed/tp1_clusters.txt
	tp2_data=$INPUTDIR/processed/tp2_clusters.txt
	source_data=$INPUTDIR/processed/source_data.tsv
	pdata=1
else
	echo "Data not processed properly"
fi


printf "\n\n\nPart 2/4:"
if [ $pdata == 1 ]; then
	Rscript scripts/cgm_collection.R --tp1 $tp1_data --tp2 $tp2_data --heights $HEIGHT
else
	echo "Not all required files found for part 2."
fi


printf "\n\n\nPart 3/4:"
if [ $pdata == 1 ]; then
	Rscript scripts/ecc_collection.R --source $source_data --strains $STRAINS --tp1 $tp1_data --tp2 $tp2_data --heights $HEIGHT --cpus 1 --trio $PARAMS
else
	echo "Not all required files found for part 3."
fi


printf "\n\n\nPart 4/4:"
if [ -f "results/ECCs.tsv" -a -f "results/CGM_strain_results.tsv" -a -f $STRAINS ]; then
  Rscript scripts/merge_data.R --ECCs results/ECCs.tsv --CGMs results/CGM_strain_results.tsv --strains $STRAINS
else
  echo "Not all required files found for part 4."
fi

printf "\nFinished process.\n"

