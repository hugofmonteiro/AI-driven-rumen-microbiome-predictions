#! /bin/bash
#
#SBATCH --mail-user=hmonteiro@ucdavis.edu   # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                     # NOTIFICATIONS OF SLURM JOB STATUS - ALL, NONE, BEGIN, END, FAIL, REQUEUE
#SBATCH -J Pred                             # JOB ID
#SBATCH -e VariableName.j%j.err             # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o VariableName.j%j.out             # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH -c 16                               # NUMBER OF PROCESSORS PER TASK
#SBATCH --mem=500Gb                         # MEMORY POOL TO ALL CORES
#SBATCH --time=07-24:00:00                  # REQUESTED WALL TIME
#SBATCH -p high2                            # PARTITION TO SUBMIT TO

# initialize conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate jup

# cd to the folder of desire, if necessary

# paste below the code to be used:

python microbiome_predictions.py

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch

