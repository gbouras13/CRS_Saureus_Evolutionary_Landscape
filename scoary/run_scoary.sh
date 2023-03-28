
### for scoary 

# assumes conda environment called panaroo exists with panaroo installed
# assumes conda environment called scoary exists with scoary installed

# diff strain
conda activate panaroo
cd Diff_strain_analysis
mkdir -p panaroo_scoary_output
panaroo -i ./T0_Gffs/*.gff -o panaroo_scoary_output -t 8 -a core --clean-mode strict  --remove-invalid-genes
conda deactivate
conda activate scoary
scoary -g panaroo_scoary_output/gene_presence_absence_roary.csv -t scoary_trait.csv

# T0 analysis

conda activate panaroo
cd T0_analysis
mkdir -p panaroo_scoary_output
panaroo -i ./Diff_Strain_Gffs/*.gff -o panaroo_scoary_output -t 8 -a core --clean-mode strict  --remove-invalid-genes
conda deactivate
conda activate scoary
scoary -g panaroo_scoary_output/gene_presence_absence_roary.csv -t scoary_trait.csv
