# SCN2A

### Requirements:
* [R](https://www.r-project.org/) with packages tidyverse, stringr, dplyr, Hmisc, memoise, reshape2, readr, logisticPCA, ROCR, corrplot, and RColorBrewer.

### Steps to Run:
* Clone the repository, modify the [config](https://github.com/helbig-lab/SCN2A/blob/master/input.yml) file.

* In the [config file](https://github.com/helbig-lab/SCN2A/blob/master/input.yml) determine the main output_dir, this is where your output files would be written to.  Default parameters have been added for most aspects of the [config file](https://github.com/helbig-lab/SCN2A/blob/master/input.yml). Processing methods and algorithms of calculation for similarity analysis can be specified by the user in the terminal after running the below script. 

* Run [R file](https://github.com/helbig-lab/SCN2A/blob/master/master_config.R), specifying the YAML config file using the --input flag .

```
~/Rscript master_config.R --input /path_to/input.yml
```

### Running the tests
There are test files available here: [Files](https://github.com/helbig-lab/SCN2A/tree/master/raw_files). Ensure that these files are linked appropriately in the [config file](https://github.com/helbig-lab/SCN2A/blob/master/input.yml) as such:

```

variant_file : raw_files/SCN2A_full.csv 

```

This provides the necessary cohort of patients with annotated HPO terms.

Note that there is an option to manually perform term propagation (see below), though default files have been provided.

## Phenotypic Similarity Analysis
Using the Human Phenotype Ontology (HPO) and a cohort of individuals annotated via HPO terms and VCF files, these scripts find phenotype-genotype correlations. This can aid in gene discovery, treatment, and a better understanding of genes' phenotypic variability. First, HPO terms' similarity scores are found for every individual pair using either the Resnik or Cube algorithm. Next, genes with potentially causitive variants in multiple individuals are extracted and the median similarity scores among each of these patients is calculated. Using permutation analysis of median similarity scores, p-values are assigned to each of these genes. A lower p-value (p < 0.05) potentially indicates a causitive gene.

### Warning
* Running the similarity analyses scripts outside of a cluster or powerful computer may take an excessively long time.  It is recommended that you run the similarity analyses in a cluster or reduce total data in your copy of the example variant file to speed up the process.

* If the ```gene_count_cube_auto.R``` file does not run, confirm that extra column was not created during initial ```cube_scn2a.csv``` file processing.  If created, the extra row can be deleted manually or the commented out line in script can be activated to remove column this automatically.

## PCA-ROC Analyses
Because there are some many potential information components that may explain the variance in phenotypic traits in individuals, these anlyses help narrow down to the most influential components in determining these phenotypic traits. These logistic analyses run in two parts.  First, the ```run_pheno_k3.R``` and  ```run_variant_k3.R``` save the two PCA models of type ```.RData```, and then the other two scripts can be run for the downstream analysis and for generating the output data frames using these models.  

## Frequency Analyses
Using the Human Phenotype Ontology (HPO) and a cohort of individuals annotated via HPO terms and VCF files, these scripts find frequency of phenotypic terms by specific groupings. For example, this analysis could be used to determine the frequency of occurrences of the HPO term 'seizures' in individuals with missense variants.  These frequencies are compared to frequencies in the remainder of each group to determine likelihood of occurence within that group. This helps determine which HPO terms may be causitively associated with a specific group. 

## Optional - Term Propagation
Although not necessary to run these scripts, as example CSVs are already provided, we've included a term propagation script within the [raw files directory](https://github.com/helbig-lab/SCN2A/tree/master/raw_files). This allows users to generate the base and propagated HPO terms on their own in order to view the process first hand. To run the propagation scripts:

* In order to create manual base and propogation files, change the yaml default parameters ``` pos_ic : raw_files/post_IC.csv ``` to ``` pos_ic :  ``` and then run the [R file](https://github.com/helbig-lab/SCN2A/blob/master/master_config.R) (see below) to generate a manual csv file.

```
~/Rscript master_config.R --input /path_to/input.yml
```

Note that this does not create the positive and negative propagation files, which are already provided in the [raw files directory](https://github.com/helbig-lab/SCN2A/tree/master/raw_files).


