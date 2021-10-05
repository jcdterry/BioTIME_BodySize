# BioTIME_BodySize

## Contents

Analysis scripts and raw data to support Terry et al. (2021) *No pervasive size trend in global community dynamics*

### Scripts

There are 5 principal scripts that generate the results in the paper from various external databases. These are .rmd files, and are accompanied with 'knitted' .html outputs that include the results within them (in the `KnittedScripts/` folder).

1. Generates assemblages from the raw BioTIME data by grouping dispersed studies into cells
2. Cleans the names assigned to records in BioTIME by cross-referencing with GBIF taxonomic backbone
3. Gathers trait data for these species from a variety of sources 
4. Conducts the analyses described in the paper

The key functions to transform the population data and calculate the slopes as described in the paper are stored in `FitThreeApproaches.R`

`ExtendedData.rmd` makes additional cosmetic changes to generate the extended data tables and figures. 

### Data Tables

Core data files that are archived elsewhere (the BioTIME database of community dynamics, the database of amniote life history traits, the fish database) are not re-hosted here, but are linked to in the scripts. Equally we do not include the raw trait data downloaded from TRY or WoRMS. 

For reproducibility of our work, we do however include the final trait tables as `bt_names_all_traits.csv` and the BioTIME metadata table for ease of reference. The derived files `bt_grid_collate_filter.csv` and `bt_grid_collate_filter_tidy.csv` are compressed into a .zip to fit on GitHub.

We include the various downloaded tables of species record names that we used in the `NameTables` folder.

`Study_Corr_Predictors.csv` details the final '$\tau$' values associates with each trait-study combination as presented in the principal figures.


## Authorship, Reuse and Licensing

Code in this repository was written by Chris Terry (Queen Mary University of London), with some parts derived from other authors as mentioned in-line. 

Feel free to use the code in this repository any way you wish, although a citation to the accompanying paper would always be welcome if relevant. Do be aware that I only attempt to link traits to species from communities with a large amount of data (as detailed in the paper), so other users may wish to repeat the process and conduct additional cleaning if they are asking different questions with a lower data requirement threshold. Equally, taxonomy is a moving target and databases are always improving. 

However, note that much of the original data used is subject to some kind of restrictions, normally the requirement to cite back to the original data compilations. Please do consult their requirements for data-reuse. The 'raw' data included in this repository is for the purpose of reproducibility and in most cases other users should return the original sources. 

