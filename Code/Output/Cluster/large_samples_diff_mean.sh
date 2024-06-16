#!/bin/bash
#SBATCH --job-name=largesamples.sh
#SBATCH --output=largesamples.out
#SBATCH --export=ALL
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20000
#SBATCH --time=10:00:00
#SBATCH --mail-user=benjamin.herzberger@student.uni-tuebingen.de
#SBATCH --mail-type=BEGIN,END,FAIL

module load math/R/4.3.3-mkl-2022.2.1-gnu-13.3 
R CMD BATCH --nosave --no-restore Large_samples_diff_mean.R
