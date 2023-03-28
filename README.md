CRS_Saureus_Evolutionary_Landscape
====================================
This repository holds the code base for Houtak et al, 'The Intra-Host Evolutionary Landscape And Pathoadaptation Of Persistent *Staphylococcus aureus* In Chronic Rhinosinusitis' (doi TBA).
------------------------------

Table of Contents
-----------
- [Assemblies](#Assemblies)
- [Chromosome Analysis Snakemake Pipeline](#Chromosome_Analysis_Snakemake_Pipeline)

- [Citation](#citation)


Assemblies
-----------

It is assumed you are using the assembled chromosome and plasmid assemblies that can be found in [PRJNA914892](https://www.ncbi.nlm.nih.gov/bioproject/914892) and also in this repository. A full list of isolates, Biosample numbers and associated metadata (particularly time-points), can be found in Supplementary Table 1. 

If you would like to recreate the assemblies, they were created with a hybrid bacterial assembly pipleine that has been formalised in a command line tool called [hybracter](https://github.com/gbouras13/hybracter). Please see the hybracter [repository](https://github.com/gbouras13/hybracter) for more details. 

Hybracter was run as follows:

```

hybracter run --input metadata.csv --output CRS_landscape_assemblies_out --threads 16

```

Where metadata.csv contains the paths to all 68 isolate long & short read FASTQ files as follows (2500000 being the lower bound for S aureus chromosome size):

C1,C1_long_read.fastq.gz,2500000,C1_short_R1.fastq.gz,C1_short_R2.fastq.gz 

These can be downloaded from the SRA. 


Chromosome Analysis Snakemake Pipeline
-----------

The section forms the bulk of the analysis conducted for the manuscript. The Snakemake pipeline can be found in teh Chromosome_Snakemake directory.

Before this, all chromosome assemblies were annotated with bakta (not included in this repository as the scripts are part of a larger project) using e.g. for C1, where C1.fasta is the chromosome assembly from hybracter:

```
bakta --db bakta_db --verbose --output C1 --prefix C1 --locus-tag C1 --threads 8 C1.fasta
```

All gff and gbk and FASTA files have been provided in this repository.

The following analyses were conducted:

1. [Snippy](https://github.com/tseemann/snippy) was run on the T1 isolate short reads vs the T0 isolate chromosome gbk for each pair to delete SNPs.
2. [Nucdiff](https://github.com/uio-cels/NucDiff) was run on the T1 isolate chromosome assembly  vs the T0 isolate chromosome assembly for each pair to detect strucutral variants.
3. [MLST](https://github.com/tseemann/mlst) was run on the T1 isolate chromosome assembly  vs the T0 isolate chromosome assembly for each pair.
4. [ISEScan](https://github.com/xiezhq/ISEScan) was run on all isolates
5. [Panaroo](https://github.com/gtonkinhill/panaroo) was run on all isolates gff files to create a pan-genome
6. [Abricate](https://github.com/tseemann/abricate) was run on all isolates to detect AMR and virulence factor genes. 
7. [Sniffles](https://github.com/fritzsedlazeck/Sniffles) was run on the T1 long reads  vs the T0 isolate chromosome assembly for each pair to detect strucutral variants.
8. [PhiSpy](https://github.com/linsalrob/PhiSpy) was run on each chromosome assemlby to predict prophages in each isolate.


To re-run these analyses, ensure you are in the Chromosome_Snakemake directory and make sure the paths to the long and short read FASTQ files (from the SRA) are changed in metadata.csv, then run:

```
snakemake -c <cores> -s runner.smk --use-conda  --conda-frontend conda  \
--config csv=metadata.csv Output=Chromosome_Analysis_Output
```



