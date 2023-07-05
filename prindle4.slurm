#!/bin/bash
# Example SLURM job script for job arrays
#
#
# This job runs a command called myprog and takes input from a set of files called 
# nnn.input, where nnn is the array task IDs
#
# set up a job array with tasks numbered 1 to 100.
#SBATCH --array=1-100
#
# give the array a single job name
#SBATCH -J prindle4
#
# Standard output is saved in this file
#SBATCH -o prindle4%A_%a.out
#
# Standard error messages are saved in this file
#SBATCH -e prindle4_%A_%a.err
#
# request one core per task
#SBATCH -c 1
#
# request a run-time of 5 mins
#SBATCH -t 2-00:00:00
#
#
# Tell SLURM if you want to be emailed when your job starts, ends, etc.
# Currently mail can only be sent to addresses @ncl.ac.uk
#SBATCH --mail-type=ALL
#
# use the task ID to locate the input file for the task.
module load MATLAB/2020a
matlab -batch "input_idx=${SLURM_ARRAY_TASK_ID}; prindle4"
