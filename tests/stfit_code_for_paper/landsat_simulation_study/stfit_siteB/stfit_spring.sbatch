#!/bin/bash
#SBATCH --nodes=1 #request one node
#SBATCH --cpus-per-task=10 #ask for 8 cpus
#SBATCH --time=10:00:00 #ask that the job be allowed to run for 2 hours.
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mingsnu@gmail.com   # Where to send mail
#SBATCH --error=job.%J.err # tell it to store the output console text to a file
#SBATCH --output=job.%J.out #tell it to store the error messages to a file
module load r #load the most recent version of R available
Rscript stfit_spring.R
