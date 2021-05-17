Metrics-CGM-ECC
================

Merged work on the two separate projects (status before initial merge
can be found at: 
* <https://github.com/vjayaman/EQProject> and 
* <https://github.com/vjayaman/ClusterGrowthMetrics>)

Testing work will be compiled and added to 
* <https://github.com/vjayaman/Testing-Metrics-CGM-ECC>

(For latest WIP - not ready for use - see *Method 3* below)

## Method 1 (Working development version):

1. Download latest (development) version from **main** branch

2. Navigate to the *Metrics-CGM-ECC/* directory in a terminal window

3. Run **environment_setup.R**, to install required R packages and checks that they can be loaded.

        .../Metrics-CGM-ECC/$ Rscript environment_setup.R

4. Open **run_process.sh** and update the input arguments to use to your data (which should be found in the *inputs/* directory)

5. Give the script *run_process.sh* executable permission

        .../Metrics-CGM-ECC/$ chmod u+x run_process.sh

6. Run the script *run_process.sh* in a terminal

        .../Metrics-CGM-ECC/$ bash run_process.sh

## Method 2 (Stable, verified version):

1. Download version [v1.0](https://github.com/vjayaman/Metrics-CGM-ECC/releases/tag/v1.0)

2. Run the following in a terminal. **Note:** no arguments need be given if your files are labelled as in the sample data (these are the default names, just replace them with your data).

    (a) **Rscript environment\_setup.R**

    This installs any R packages required but not found in the environment, and checks that they can be loaded without issue.


    (b) **Rscript prepare\_inputs.R** -i *\<input directory\>* -m *\<metadata filename\>* -a *\<TP1 cluster assignments - unprocessed\>* -b *\<TP2 cluster assignments - unprocessed\>*

    To run, it needs 3 files in the *inputs/* directory: a cluster assignments file for time point 1, a cluster assignments file for time point 2, and a strain metadata file, all tab-delimited. 
    
    It also writes a placeholder file for the source distances file to the newly created *inputs/processed/* directory. 
    
    To use files with the same names and location as the sample data, simply run:

        $ RScript prepare_inputs.R
    
    or to manually feed input files:

        $ Rscript prepare_inputs.R -i inputs -m strain_info.txt -a tp1_clusters_init.txt -b tp2_clusters_init.txt

    
    (c) **Rscript CGM/datacollection.R** -a *\<tp1 cluster assignments\>* -b *\<tp2 cluster assignments\>* -x *\<comma-delimited string of heights, or a single height\>*

    If you do not change the filenames in *inputs/processed/*, you can run this like:

        $ Rscript CGM/datacollection.R
    
    or to manually feed input files: 

        $ Rscript CGM/datacollection.R -a inputs/processed/tp1_clusters.txt -b inputs/processed/tp2_clusters.txt -x 0

    This collects the cluster growth metrics and saves them to newly created *results/* directory.

    (d) **Rscript ECC/cov\_epi.R** -a *\<source data\>* -b *\<strain data\>* -c *\<TP1 cluster assignments\>* -d *\<TP2 cluster assignments\>* -x *\<comma-delimited string of heights\>* -t *\<hyphen-delimited numeric codes for source, temporal, geographic coefficients\>*

    Note that the source data file does not need to be provided, the script
will use the default placeholder if this is left blank. This saves the
resulting ECCs to the *results/* directory. The following example
collects ECCs for height 0 with coefficients of source = 0, temporal =
1, and geo = 0:

        $ Rscript ECC/cov_epi.R -a inputs/processed/source_data.tsv -b inputs/strain_info.txt -c inputs/processed/tp1_clusters.txt -d inputs/processed/tp2_clusters.txt -x "0" -t "010"

    Another example, but using the placeholder source data (since the
coefficient is 0 anyway) and for two sets of ECC parameters (source = 0,
temporal = 0, and geo = 1) and (source = 0, temporal = 1, and geo = 0):

        $ Rscript ECC/cov_epi.R -b inputs/strain_info.txt -c inputs/processed/tp1_clusters.txt -d inputs/processed/tp2_clusters.txt -x "0" -t "001-010"

    To merge the ECCs and the CGM results into a strain file and a cluster
file:

    (e) **Rscript merge\_data.R -e** *\<ECC results\>* **-c** *\<CGM results\>*
**-s** *\<Strain metadata\>*

    As with the other files, the input arguments for this can be left blank to use the default file names.
    
        $ Rscript merge_data.R -e results/ECCs.tsv -c results/CGM_strain_results.txt -s inputs/strain_info.txt

## Method 3 (Latest version, not ready for use):

* The branch `dist_mat_memory` contains the refactored and latest version, which deals with issues relating to large input datasets (> 35000 samples, for example).

* Before, a distance matrix of the required size couldn't be held in memory ("R cannot allocate vector of size ...").

* This has been dealt with, but there are still a few bugs/issues that need to be addressed:
	* Type handling has been done for the CGM part, but ECC inheritance needs to be re-done (more efficiently)
	* the average distance columns (from the ECC results) is much too slow - and memory intensive for sizeable datasets --> need a better method than going cluster by cluster and merging each time
	* need to bring back in the original clusters given (if they were not numbers in initial input) as another column or two (since they are represented by numbers during analysis)
	* need to use quosures so the deltaECC columns can be made without referring to the specific parameters
	* and any modifications that I overlooked from before (including the more detailed input parameter file)

* I can stitch the new method back into the pipeline then and then finish fine-tuning the output tables changes mentioned in the last meeting.

* This will be merged back into the main branch, and after someone else has verified that everything still runs for them as expected, this will be tagged as *v2.0*.

