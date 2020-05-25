#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=6000M
#SBATCH --account=rrg-cgreenwo
#SBATCH --mail-user=kaiqiong.zhao@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END


module load r
echo Started: 
date
Rscript "/project/6002088/zhaokq/Simulation/Exp_Null_Z_BSmooth/Simu_1000_Samp_100_BSmooth.R" 
echo Ended
date


