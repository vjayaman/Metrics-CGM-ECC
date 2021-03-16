# Metrics-CGM-ECC

Merged work on the two separate projects (status before initial merge can be found at: https://github.com/vjayaman/EQProject and https://github.com/vjayaman/ClusterGrowthMetrics)

To collect the results, run the following in a terminal:

  (1) Rscript environment_setup.R

  (2) Rscript prepare_inputs.R

  (3) Rscript CGM/datacollection.R -a inputs/processed/tp1_clusters.txt -b inputs/processed/tp2_clusters.txt -x "0"

  (4) Rscript ECC/cov_epi.R -a inputs/processed/source_data.tsv -b inputs/strain_data.tsv -c inputs/processed/tp1_clusters.txt
-d inputs/processed/tp2_clusters.txt -x "0" -t "010-001"

  (5) Rscript merge_data.R -e results/ECCs.tsv -c results/CGM_strain_results.txt -s inputs/strain_data.tsv

