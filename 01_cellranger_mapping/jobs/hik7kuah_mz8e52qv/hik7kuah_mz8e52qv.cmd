#!/usr/bin/env bash

#SBATCH --job-name="hik7kuah_mz8e52qv"
#SBATCH --workdir=.

#SBATCH --error=./logs/slurm_%x_%J.err
#SBATCH --output=./logs/slurm_%x_%J.out

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --partition=genB,main

#SBATCH --mail-type=END       
#SBATCH --mail-user=laura.jimenez@cnag.crg.eu

echo [`date "+%Y-%m-%d %T"`] started job on $HOSTNAME

ulimit -n 16000
export HDF5_USE_FILE_LOCKING="FALSE"

/scratch/groups/singlecell/software/cellranger/6.0.1/cellranger multi --id hik7kuah_mz8e52qv \
    --csv /scratch/devel/ljimenez/projects/POST/POST_12/01_cellranger_mapping/jobs/hik7kuah_mz8e52qv/hik7kuah_mz8e52qv_config.csv \
    --localcores 2 \
    --jobmode /scratch/groups/singlecell/software/cellranger/6.0.1/external/martian/jobmanagers/slurm.template

echo [`date "+%Y-%m-%d %T"`] finished job
