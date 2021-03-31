#!/usr/bin/env bash

INPUTDIR=inputs/
STRAINS=strain_info.txt
TP1_CLUSTERS=tp1_clusters_init.txt
TP2_CLUSTERS=tp2_clusters_init.txt

printf "Part 1/4:\n"
Rscript prepare_inputs.R --inputdir $INPUTDIR --metadata $STRAINS --tp1 tp1_clusters_init.txt --tp2 tp2_clusters_init.txt
printf "\n\n\n"

printf "Part 2/4:\n"
Rscript CGM/datacollection.R --tp1 $INPUTDIR/processed/tp1_clusters.txt --tp2 $INPUTDIR/processed/tp2_clusters.txt --heights 0
printf "\n\n\n"

printf "Part 3/4:\n"
Rscript ECC/cov_epi.R --source $INPUTDIR/processed/source_data.tsv --strains $INPUTDIR/strain_info.txt --tp1 $INPUTDIR/processed/tp1_clusters.txt --tp2 $INPUTDIR/processed/tp2_clusters.txt --heights 0 --cpus 1 --trio "010-001"
printf "\n\n\n"

printf "Part 4/4:\n"
Rscript merge_data.R --ECCs results/ECCs.tsv --CGMs results/CGM_strain_results.tsv --strains $INPUTDIR/$STRAINS
printf "\n\n\n"

echo 'Finished process.'
