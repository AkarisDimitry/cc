#!/bin/bash
#
#SBATCH --nodes=1    
#SBATCH -J {jobname}
#SBATCH --time={time_cap}
#SBATCH --tasks-per-node={ntasks}
#SBATCH -p {partition}
#

# load modules
module load 2023
module load VASP6/6.4.2-foss-2023a-VASPsol-VTST

# launch VASP
srun vasp_std >>out 2>>err