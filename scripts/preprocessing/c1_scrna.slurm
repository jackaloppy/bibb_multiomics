#!/bin/bash

#SBATCH --job-name=Bibb-c1-scrna-hu
#SBATCH --nodes=1  
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=5gb            
#SBATCH --time=2:00:00   
#SBATCH --account=xiaosun
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_xiaosun
 
module load cellranger/7.0.0

cd /xdisk/yiwenliu/ranwei/bibb
   
cellranger count --id=C1-scRNA-Human \
   --fastqs=/xdisk/xiaosun/ranwei/bibb/scrna \
   --sample=C1 \
   --transcriptome=/groups/xiaosun/ranwei/ra_seq_proj/run_cellranger_count/refdata-gex-GRCh38-2020-A \
   --chemistry=ARC-v1