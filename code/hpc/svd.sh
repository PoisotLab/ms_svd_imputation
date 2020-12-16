#! /bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=virus-svd
#SBATCH --array=1-912%50
#SBATCH --output=%x-%a.out
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2300M

module load StdEnv/2020 julia/1.5.2

julia --project -t 39 main.jl
