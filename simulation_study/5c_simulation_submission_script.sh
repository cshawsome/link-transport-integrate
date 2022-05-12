#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o /u/home/c/cshaw343/joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=1:00:00,h_data=4G #requests 1 hours, 2GB of data (per core) 
#$ -pe shared 1 #requests 1 core
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# TO TEST MATTER RUN ONLY TWICE:
##$ -t 1-2:1
# FOR THE FULL RUN USE INSTEAD:
#$ -t 1-45:1
## 

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.1.0 #loads R/4.1.0 for use 
export OMP_NUM_THREADS=1 #uses max 1 thread (needs to match -pe shared)
## 
# run R code

#The jobs are combinations of
#
#distribution = c("normal", "lognormal", "bathtub")
#sample_size = c(500, 1000, 2000, 4000, 8000)
#prior_prop = c("ADAMS", "unimpaired", "dementia")
#
#the task ID indicates the row of the simulation scenario .csv file 

echo "======"
echo SGE_TASK_ID=$SGE_TASK_ID      
R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID num_replicates=10"  /u/home/c/cshaw343/RScripts/simulation_run.R /u/home/c/cshaw343/output/output.$JOB_ID.$SGE_TASK_ID
echo R CMD BATCH --no-save --no-restore \'--args scenario_num=$SGE_TASK_ID num_replicates=10 \'  /u/home/c/cshaw343/RScripts/simulation_run.R /u/home/c/cshaw343/output/output.$JOB_ID.$SGE_TASK_ID

