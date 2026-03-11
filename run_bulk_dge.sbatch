#!/bin/bash
#SBATCH -J bulk_dge
#SBATCH --time=0-4:00:00
#SBATCH --partition=b64_1024
#SBATCH --exclude=b[048-049]
#SBATCH -N 1 		#set the desired number of nodes
#SBATCH --mem=100G 	#memory amount has to be specified otherwise the whole amount will be dedicated to this job and the other one won't be able to enter the node
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=me@lab.org
#SBATCH --mail-type=END,FAIL

source /uochb/soft/b/spack/20221211-git/share/spack/setup-env.sh
echo "Loaded singularity from spack"

source /uochb/soft/generic/miniconda/202310/bin/activate # get miniconda
conda init
echo "Conda activated"

conda activate /home/faltejskova/.conda/envs/bulkdge

python some/path/bulk_rnaseq_DGE.py --jobname bulk_dge --user-email me@lab.org --count-matrix counts.tsv --sample-file samples.tsv --contrasts contrasts.tsv --gtf-file genes.gtf --padj-alpha 0.05 --fold-change-threshold 1 --output-directory-path bulk_dge_analysis_output --minimal-gene-count 15 --interesting-genes-file interesting

