#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=6000M
#SBATCH --account=rrg-cgreenwo
#SBATCH --mail-user=kaiqiong.zhao@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --array=1-14

echo Started: 
date
echo $SLURM_ARRAY_TASK_ID
Rscript "/project/6002088/zhaokq/Simulation/Run-BS-SMSC-No-Error/run_100.R" $SLURM_ARRAY_TASK_ID
echo Ended
date
