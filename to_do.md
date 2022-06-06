
## To do

### Last checked: May 16, 2022

#### Done

* change date formatting so Excel doesn’t convert weeks (force as characters) – uses YearWeek instead of just Week
* output a table of || TP2 cluster | Pango lineage | Size ||
    * check this is easy to extract for user
* want to minimize user text entry
* should be able to get top 10 clusters by size and by rate (ones with fastest growth rate)
* clustering before the heatmap 
    * should be clearly clustered
* pull out the heatmap color palette (should be user-adjustable)
* cluster by (geo, temp, temp-geo) – user choice
* the whisker plot run out of memory if running all top 20 clusters when there are two cluster sizes bigger than 9000 (comment out)
* heatmaps – original cluster name 
* heatmap adjustment – for largest growth clusters (what is minimum size to consider? note: min value of 2, will not consider singletons)
* weekly ECC not consistent with monthly ECC at turn of year
    * new report missing the first week of 2021, use YearWeek instead of Week
* transform in heatmap section should be off as default
    * add as option in form_inputs
* novel_only ECC columns not found
* add options for YearWeek labeling (years x, y, where x is the year before y) in the case of overlapping between last week of x and first week of y - suppose last week of x is week x-52: 
    * first week of y is labeled x-52
    * first week of y is labeled y-01
* paired histograms – should be able to shade these
    * see scripts/Misc/plotting_functions.R and the frequencyPlots() function for how it's plotted
* documents describing the analysis – figure out which components need this (of my part)
* make sure the notes for clusters with figures are clear
    * see *results/<interval type>/cluster_labs.txt* after running script 7 to get details on which clusters have figures
* top growth clusters not the same as GZ’s results
    * fixed --> see *all_reports/Epiquant_Report_05-18-2022.Rmd* for the check on cluster sizes. All of the report now works with the latest outputs from the metrics generation part, and the document formatting now has more spacing, etc. for better readability
* line 919 – pairwise geo and temp distance for top 20 clusters need to calculate the actual physical pairwise distance, then we have to change the distance calculation in step 4
    * they are calculating actual distances in the metrics generation process (but need to articulate this better in the docs)

#### Not done

* better average_dists implementation (faster, for larger files) (use C/python)
* better ECC collections (faster, for larger files, consider parallelization) (use C)
* whisker plots – too much data?
* heatmap generation - write up in a jupyter notebook / Rmarkdown file - so it's easy for others to verify steps
    * partially done
* initial setup should prompt user to remove/rename *logs/*, *intermediate_dists/*, and *results/* if running a new set of data (so overwriting results for a different data set does not happen)

