#!/bin/bash

#SBATCH --job-name=Bibb-c1-arc
#SBATCH --nodes=1  
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=5gb            
#SBATCH --time=12:00:00   
#SBATCH --account=xiaosun
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_xiaosun
 
export PATH=/groups/xiaosun/ranwei/ra_seq_proj/run_cellranger_count/cellranger-arc-2.0.2:$PATH

cd /xdisk/xiaosun/ranwei/bibb

cellranger-arc count --id=C1 \
   --reference=/groups/xiaosun/ranwei/ra_seq_proj/run_cellranger_count/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
   --libraries=/home/u27/ranwei/ra/bibb_code/libraries.csv