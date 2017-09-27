#!/bin/bash
#SBATCH --job-name sync2box
#SBATCH --account=IC1066_boppana            ##Billing account
#SBATCH --ntasks=1                          ##Number of PROCESSES
#SBATCH --cpus-per-task=1                   ##Number of CORES per PROCESSS
#SBATCH --mem-per-cpu=2000                  ##Memory specified for each core used (in MB) (no cores, use --mem=)
#SBATCH -t 0-12:00:00 --partition=short     ##Runtime in D-HH:MM:SS
#                                           ## express(2h), short(12h), medium(2d2h), long(6d6h), interactive(2h)
#SBATCH --share
#
#SBATCH --mail-user=%u@uab.edu              ## assume blazerID
#SBATCH --mail-type=ALL                     ## BEGIN, END, ERROR, ALL
#
#SBATCH --error=logs/%j.%N.err.txt               ##File to which STDERR will be written
#SBATCH --output=logs/%j.%N.out.txt              ## File to which STDOUT will be written
srun ~/uab_ngs/box/sync2box
