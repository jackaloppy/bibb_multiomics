#!/bin/bash

#SBATCH --job-name=Bibb-c1-scatac
#SBATCH --nodes=1  
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=5gb            
#SBATCH --time=6:00:00   
#SBATCH --account=xiaosun
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_xiaosun
#SBATCH --mail-type=END
 
export PATH=/groups/xiaosun/ranwei/ra_seq_proj/run_cellranger_count/cellranger-atac-2.1.0:$PATH

cd /xdisk/yiwenliu/ranwei/bibb
   
cellranger-atac count --id=C1-scATAC \
   --fastqs=/xdisk/xiaosun/ranwei/bibb/scatac/C1-ATAC \
   --sample=C1 \
   --reference=/groups/xiaosun/ranwei/ra_seq_proj/run_cellranger_count/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
   --chemistry=ARC-v1