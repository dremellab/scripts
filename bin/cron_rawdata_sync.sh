#!/bin/bash
#SBATCH --job-name=rawdata_sync           # Default job name
#SBATCH --time=12:00:00                   # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=1                 # Cores per task
#SBATCH --mem=40G                         # Memory
#SBATCH --partition=standard              # Partition name
#SBATCH --account=dremel_lab    	  # account

set -exo pipefail

rsync -avh --update /standard/dremel_lab/rawdata/ /scratch/$USER/rawdata/
rsync -avh --update /scratch/$USER/rawdata/ /standard/dremel_lab/rawdata/

# resubmit for tomorrow
sbatch --begin=tomorrow "$0"
