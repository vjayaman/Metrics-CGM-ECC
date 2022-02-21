WIP - under development









Imported to GitLab - cscscience setup. Will synchronize updates here with the GitLab repository with each major change. 

Metrics-CGM-ECC
================

## Method 1 (Latest version, WIP):

Updated: Feb 21, 2022

1. Download latest version (will include more details when it is refined enough to be tagged and released)

2. Move the input files you want to use (a strain metadata file and one cluster assignment file - last time point, after all strains have been added) into *Metrics-CGM-ECC/inputs/*, then rename to match the examples in the *inputs/* directory.
 
3. Open **Metrics-CGM-ECC/form.html** in a browser, then fill out the form (arguments should correspond to your input file data)
	* Download the resulting text file, and move it to the *inputs/* directory

4. The next section can be run in one of three ways

### Method 1a

5. Navigate to the *Metrics-CGM-ECC/* directory in a terminal window. 

6. Run **environment_setup.rscript**, to install the required R packages and check that they can be loaded.

        .../Metrics-CGM-ECC/$ Rscript environment_setup.rscript

7. Give the script *collect_metrics.rscript* executable permission

        .../Metrics-CGM-ECC/$ chmod u+x collect_metrics.rscript

8. Run the script *collect_metrics.rscript* in a terminal window 

        .../Metrics-CGM-ECC/$ Rscript collect_metrics.rscript
	
### Method 1b

5. The scripts can be run individually in RStudio. Open up RStudio and use *setwd()* to navigate to the Metrics-CGM-ECC/ directory.

6. Use *source('environment_setup.rscript')* in the console or use Code > Source after opening *environment_setup.rscript*

7. Use the same idea as in 6. for all numbered files in Metrics-CGM-ECC/scripts/, from 1\_prepare\_inputs.rscript to 7b\_topY\_heatmaps.rscript

### Method 1c

5. Open up the Metrics-CGM-ECC/ directory using something like a file manager. 

6. Double-click on *environment_setup.rscript* - you may be prompted by a security window to give the script run permissions. Please do so (you can open the file using notepad or RStudio, for example, to verify that nothing risky is being done)

7. After the new directory structure is settled (you can see a logs/ folder, for example), double-click on *collect_metrics.rscript*

The CGM data, ECC metrics, and merged files (for strains and clusters) will be saved to a newly created *results/* directory.


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
	

## Notes:

Merged work on the two separate projects (status before initial merge
can be found at: 
* <https://github.com/vjayaman/EQProject> and 
* <https://github.com/vjayaman/ClusterGrowthMetrics>)

Testing work will be compiled and added to 
* <https://github.com/vjayaman/Testing-Metrics-CGM-ECC>
