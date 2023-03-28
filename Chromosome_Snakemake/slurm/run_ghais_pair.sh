#!/bin/bash -l

#SBATCH --job-name=ghais_pair
#SBATCH --mail-user=george.bouras@adelaide.edu.au
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --err="ghais_pair.err"
#SBATCH --output="ghais_pair.out"

# Resources allocation request parameters
#SBATCH -p batch
#SBATCH -N 1                                                    # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 32                                                    # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=12:00:00                                         # time allocation, which has the format (D-HH:MM), here set to 1 hou                                           # generic resource required (here requires 1 GPUs)
#SBATCH --mem=75GB                                              # specify memory required per node


# run from Bacteria_Multiplex


PROF_DIR="/hpcfs/users/a1667917/snakemake_slurm_profile"

# move up a directory

cd ..

module load Anaconda3/2020.07
conda activate snakemake_clean_env

# snakemake -c 1 -s runner_ghais.smk --use-conda  --conda-frontend conda --conda-create-envs-only  \
# --config csv=ghais_metadata.csv Output=../Paper_1_Snakemake_Output 

snakemake -c 32 -s runner_ghais.smk --use-conda  --conda-frontend conda  \
--config csv=ghais_metadata.csv Output=../Paper_1_Snakemake_Output 




conda deactivate
