Metrics-CGM-ECC
================

Latest verified version: v1.0

Merged work on the two separate projects (status before initial merge
can be found at: <https://github.com/vjayaman/EQProject> and
<https://github.com/vjayaman/ClusterGrowthMetrics>)

To collect the results, run the following in a terminal:

**Rscript environment\_setup.R**

This installs any R packages required but not found in the environment,
and checks that they can be loaded without issue. It also writes a
placeholder source distances file to the *inputs/processed/* directory,
which can be overwritten if you have source data you want to use
instead.

If running this process on the sample data, run

**Rscript prepare\_inputs.R**

otherwise move your (*correctly formatted*) input files into the
*inputs/processed/* directory. To run the following, 4 tab-delimited
files are needed: a cluster assignments file for time point 1, a cluster
assignments file for time point 2, a strain metadata file, and a source
distances file (will use the placeholder if one not provided).

Example:
`Rscript prepare_inputs.R -i inputs -m strain_info.txt -a tp1_clusters_init.txt -b tp2_clusters_init.txt`

If the sample data is used, and `prepare_inputs.R` is run, you get the
files used in the examples below.

**Rscript CGM/datacollection.R -a** *\<tp1 cluster assignments\>* **-b**
*\<tp2 cluster assignments\>* **-x** *\<comma-delimited string of
heights\>*

Example:  
`Rscript CGM/datacollection.R -a inputs/processed/tp1_clusters.txt -b
inputs/processed/tp2_clusters.txt -x "0"`

This collects the cluster growth metrics and saves them to newly created
*results/* directory.

**Rscript ECC/cov\_epi.R -a** *\<source data\>* **-b** *\<strain data\>*
**-c** *\<TP1 cluster assignments\>* **-d** *\<TP2 cluster
assignments\>* **-x** *\<comma-delimited string of heights\>* **-t**
*\<hyphen-delimited numeric codes for source, temporal, geographic
coefficients\>*

Note that the source data file does not need to be provided, the script
will use the default placeholder if this is left blank. This saves the
resulting ECCs to the *results/* directory. The following example
collects ECCs for height 0 with coefficients of source = 0, temporal =
1, and geo = 0:

Example:  
`Rscript ECC/cov_epi.R -a inputs/processed/source_data.tsv -b
inputs/strain_info.txt -c inputs/processed/tp1_clusters.txt -d
inputs/processed/tp2_clusters.txt -x "0" -t "010"`

Another example, but using the placeholder source data (since the
coefficient is 0 anyway) and for two sets of ECC parameters (source = 0,
temporal = 0, and geo = 1) and (source = 0, temporal = 1, and geo = 0):

Example: `Rscript ECC/cov_epi.R -b inputs/strain_info.txt -c
inputs/processed/tp1_clusters.txt -d inputs/processed/tp2_clusters.txt
-x "0" -t "001-010"`

To merge the ECCs and the CGM results into a strain file and a cluster
file:

**Rscript merge\_data.R -e** *\<ECC results\>* **-c** *\<CGM results\>*
**-s** *\<Strain metadata\>*

Example: `Rscript merge_data.R -e results/ECCs.tsv -c
results/CGM_strain_results.txt -s inputs/strain_info.txt`
