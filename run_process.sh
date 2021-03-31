#!/usr/bin/env bash

# USER-SPECIFIED INPUTS, change these if necessary before running:
INPUTDIR=inputs
STRAINS=$INPUTDIR/strain_info.txt
TP1_CLUSTERS=$INPUTDIR/tp1_clusters_init.txt
TP2_CLUSTERS=$INPUTDIR/tp2_clusters_init.txt

# Do not change these variables unless necessary:
pdata=0
SOURCE=use_placeholder


printf "\n\nPart 1/4:"
if [ -f $STRAINS -a -f $TP1_CLUSTERS -a -f $TP2_CLUSTERS ]; then
	Rscript prepare_inputs.R --inputdir $INPUTDIR --metadata $STRAINS --tp1 $TP1_CLUSTERS --tp2 $TP2_CLUSTERS --source $SOURCE
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
	Rscript CGM/datacollection.R --tp1 $tp1_data --tp2 $tp2_data --heights 0
else
	echo "Not all required files found for part 2."
fi


printf "\n\n\nPart 3/4:"
if [ $pdata == 1 ]; then
	Rscript ECC/cov_epi.R --source $source_data --strains $STRAINS --tp1 $tp1_data --tp2 $tp2_data --heights 0 --cpus 1 --trio "010-001"
else
	echo "Not all required files found for part 3."
fi


printf "\n\n\nPart 4/4:"
if [ -f "results/ECCs.tsv" -a -f "results/CGM_strain_results.tsv" -a -f $STRAINS ]; then
  Rscript merge_data.R --ECCs results/ECCs.tsv --CGMs results/CGM_strain_results.tsv --strains $STRAINS
else
  echo "Not all required files found for part 4."
fi


printf "\n\n\nFinished process."
