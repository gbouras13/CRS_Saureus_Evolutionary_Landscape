# Plasmid_Analysis_Pipeline
Pipeline to Analyse Bacterial Plasmid Sequences

This pipeline calculates the mash score between all plasmids in the dataset.

https://mash.readthedocs.io/en/latest/index.html

The output is a distance from 0 to 1, where 0 is identical, and 1 is completely different.

The pipeline also annotates the plasmids using bakta and creates a pangenome using panaroo.

### Usage

snakemake -c 1 -s plasmid_runner.smk --use-conda   \
--config Input=../PLASMID_FASTAS Output=../Plasmid_Snakemake_Out
