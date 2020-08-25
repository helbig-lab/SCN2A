# SCN2A

### Requirements:
* [R](https://www.r-project.org/) with packages tidyverse, stringr, dplyr, readr and memoise.

### Steps to Run:
* Clone the repository, modify the [config](https://github.com/helbig-lab/SCN2A/blob/master/input.yml) file.

* In the [config file](https://github.com/helbig-lab/SCN2A/blob/master/input.yml) mention the the field output_dir, this is where your output files would be written to.  Processing methods and algorithms of calculation for similarity analysis can be specified by the user in the terminal after running the below script. 

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

## Phenotypic Similarity Analysis
Using the Human Phenotype Ontology (HPO) and a cohort of patients annotated via HPO terms and VCF files, these scripts find phenotype-genotype correlations. This can aid in gene discovery, treatment, and a better understanding of genes' phenotypic variability. First, HPO terms similarity scores are found for every patient pair using either the Resnik or Cube algorithm. Next, genes with potentially causitive variants in multiple patients are extracted and the median similarity scores among each of these patients is calculated. Using permutation analysis of median similarity scores, p-values are assigned to each of these genes. A lower p-value potentially indicates a causitive gene.

### Warning
* Running this script outside of a cluster or powerful computer may take an excessively long time.  It is recommended that you run the similarity analyses in a cluster or remove data from your copy of the example variant file to speed up the process.

* If the ```gene_count_cube_auto.R``` file does not run, confirm that extra row was not created during initial ```cube_scn2a.csv``` file processing.  If created, the extra row can be deleted manually or commented out line in script can be activated.


## Frequency Analyses
Using the Human Phenotype Ontology (HPO) and a cohort of patients annotated via HPO terms and VCF files, these scripts find frequency of phenotypic terms by specific groupings. For example, this analysis could be used to determine the frequency of occurrences of the term 'seizures' in patients with missense variants.  These frequencies are compared to frequencies in the remainder of each group to determine likelihood of occurence within that group. This helps determine which HPO terms may be causitively associated with a specific group. 
