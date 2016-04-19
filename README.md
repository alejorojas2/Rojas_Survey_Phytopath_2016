Oomycete species associated with soybean seedlings in North America
=======

This repository contains data and scripts used to analyze data for the manuscripts:

1. Oomycete species associated with soybean seedlings in North America - Part I: identification and pathogenicity characterization
2. Oomycete species associated with soybean seedlings in North America - Part II: diversity and ecology in raltion to environmental and edaphic factors

Data, code and results are contained in this project, but there are independent [Rmarkdown] documents for each analysis.

Overview
--------

    project
    |- README          # the top level description of content
    |
    |- data            		# raw and primary data, are not changed once created 
    |  |- raw/         		# raw data, will not be altered
    |  |- clean/       		# cleaned data, will not be altered once created
    |  |- input_files_mothur	#Input files used for mothur analyses
    |
    |- code/           # any programmatic code
    |  |- abundance_analysis.R #Analyses part I    |  |- seed_rot_analysis.R #Analyses part I    |  |- seedling_analysis.R #Analyses part I    |  |- soy_sample_map.R #Analyses part I    |  |- tree_traits.R #Analyses part I
    |  |- mothur_field.batch # Batch file used in mothur   
    |
    |- results         # all output from workflows and analyses
    |  |- figures/     # graphs, likely designated for manuscript figures
    |  |- Oomycetes_part-I_analysis.Rmd    |  |- Seed_rot_data_analysis.Rmd    |  |- Seedling_data_analysis.Rmd
    |  |- Oomycetes_part-II_analysis.Rmd
    |  |- Ecology_analyses_oomycetes.Rmd
    |
    |- scratch/        # temporary files that can be safely deleted or lost
    |
    |- study.Rmd       # executable Rmarkdown for this study, if applicable
    |- study.Rproj     # RStudio project for this study, if applicable


Datasets
----------
For both papers, the sequences can be obtained with the accession numbers GenBank KU208091-KU211502 at [GenBank](http://www.ncbi.nlm.nih.gov/nuccore/).

###Part I
* __Isolates_11-12_final.csv__: File containing identification of isolates using local BLAST database based on [Robideau et al. 2011](http://onlinelibrary.wiley.com/doi/10.1111/j.1755-0998.2011.03041.x/abstract)
* __seed_rot.csv__: File containing seed rot ratings of soybean cv. 'Sloan' challenged with oomycete species.
* __dry_weight_allsets.csv__: File containing dry weights for soybean cv. 'Sloan' seedlings challenged with oomycete species.
* __root_measurements_final.csv__: File containing root area, root length for soybean cv. 'Sloan' seedlings challenged with oomycete species.

###Part II
* __Seqs_11_12.OTU.biom__: BIOM file containing OTU data
* __env_metadata.csv__: Environmental and edaphic data by field
* __OTU.phylotype.txt__: Phylotype data for analyses
* __Taxa.phylotype.txt__ Taxonomy assignments for phylotype

###Miscelaneous datasets
* __CA_GIS_data.txt__: Canada locations sampled in this study
* __US_GIS_data.txt__: U.S. locations sampled in this study
* __soybean_data.csv__: Soybean planted are in the U.S. in 2012 (NASS USDA).
* __Ontario_soybean.txt__: Soybean production in Ontario, CA.


Results
----------------------
Rmarkdown documents containing figures and analyses for Oomycete species associated with soybean seedlings in the U.S.

* [Part I - Identification and pathogenicity characterization](results/Oomycetes_part-I_analysis.md)

* [Part II - Diversity and ecology in relation to environmental and edaphic factors](results/Oomycetes_part-II_analysis.md)



Miscelaneous
----------------
The initial file and directory structure of this project was developed by a group of participants in the Reproducible Science Curriculum Workshop, held at [NESCent] in December 2014. 

To access this repository template use [rr-init repository](https://github.com/Reproducible-Science-Curriculum/rr-init)

[NESCent]: http://nescent.org
[Rmarkdown]: http://rmarkdown.rstudio.com/
