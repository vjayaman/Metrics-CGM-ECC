#!/usr/bin/env bash

# use sed -i 's/\r//g' practice.sh to remove Windows carriage returns
start_time=$SECONDS

# Strain metadata file, with "Strain" identifier column
rawstrains=inputs/strain_info.txt

# Unprocessed cluster assignments, with "Strain" identifier column
tp1clusters=inputs/tp1_clusters_init.txt
tp2clusters=inputs/tp2_clusters_init.txt
form=inputs/form_inputs.txt

# DO NOT CHANGE pdata, it's a flag used to check for success
pdata=0



printf "\n\nPart 0/6:"
Rscript environment_setup.R
step0=$(($SECONDS - $start_time))



printf "\n\nPart 1/6:"
if [ -f $rawstrains -a -f $tp1clusters -a -f $tp2clusters -a -f $form ]; then
	Rscript scripts/1_prepare_inputs.R
else
	echo "Not all required files found for part 1."
fi

if [ -f inputs/processed/tp1_clusters.txt -a -f inputs/processed/tp2_clusters.txt ]; then
	tp1_data=inputs/processed/tp1_clusters.txt
	tp2_data=inputs/processed/tp2_clusters.txt
	strains=inputs/processed/strain_info.txt
	pdata=1
else
	echo "Data not processed properly"
fi
step1=$(($SECONDS - $step0))



printf "\n\n\nPart 2/6:"
if [ $pdata == 1 ]; then
	Rscript scripts/2_cgm_collection.R
else
	echo "Not all required files found for part 2."
fi
step2=$(($SECONDS - $step1))



printf "\n\n\nPart 3/6:"
if [ $pdata == 1 ]; then
	Rscript scripts/3_dist_matrices.R
else
	echo "Not all required files found for part 3."
fi
step3=$(($SECONDS - $step2))



printf "\n\n\nPart 4/6:"
if [ $pdata == 1 ]; then
	Rscript scripts/4_ecc_collection.R -m $strains -a $tp1_data -b $tp2_data -d $form
else
	echo "Not all required files found for part 4."
fi
step4=$(($SECONDS - $step3))


# # printf "\n\n\nPart 5/6:"
# # if [ -f "results/ECCs.tsv" -a -f "results/CGM_strain_results.tsv" -a -f $strains ]; then
# #   Rscript scripts/5_combine_results.R -m $strains -a inputs/processed/allTP1.Rds -b inputs/processed/allTP2.Rds -c results/CGM_strain_results.tsv -d $form -e results/ECCs.tsv
# # else
# #   echo "Not all required files found for part 5."
# # fi
# # 
# # printf "\n\n\nPart 6/6:"
# # if [ -f "results/CGM_strain_results.tsv" -a -f "intermediate_data/dist_extremes.Rds" -a -f $strains ]; then
# # 	Rscript scripts/6_heatmap_dists.R -m $strains -a $tp1_data -b $tp2_data -c results/CGM_strain_results.tsv -d $form
# # else
# # 	echo "Not all required files found for part 6."
# # fi
# 

echo ""
echo "Part 0 took $(($step0 / 60)) minutes and $(($step0 % 60)) seconds."
echo "Part 1 took $(($step1 / 60)) minutes and $(($step1 % 60)) seconds."
echo "Part 2 took $(($step2 / 60)) minutes and $(($step2 % 60)) seconds."
echo "Part 3 took $(($step3 / 60)) minutes and $(($step3 % 60)) seconds."
echo "Part 4 took $(($step4 / 60)) minutes and $(($step4 % 60)) seconds."

printf "\nFinished process.\n"