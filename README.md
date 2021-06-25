# PRECISION clonality paper

R package and analysis scripts used in the PRECISION clonality paper. 

Directory structure
--------------------

* `breakclone` - R package to assess clonality of tumour pairs.
* `analysis` - Scripts used in the analysis.

Analysis
--------

This directory contains the following scripts used in the analysis:

* `plotSegs.R` R function which takes as an input a QDNAseq object and plots the copy number profile of an specific sample.
* `make_oncoPrint.R` takes as an input table of mutations and outputs statistical analysis of number of mutations for all genes in the panel and plots oncoPrint.
