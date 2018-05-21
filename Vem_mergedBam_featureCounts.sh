#!/bin/bash

#SBATCH -p physical
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=exon_count
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err 



/home/afarre/.local/necklace-necklace_v0.9/tools/bin/featureCounts -T 32 --primary -p -t exon -g gene_id \
 -a /home/afarre/Vem_necklace/mapped_reads/stringtie_merged_out.bam/stringtie_merged_bam_blocks.gtf \
  -o /home/afarre/Vem_necklace/mapped_reads/stringtie_merged_out.bam/gene.counts \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030152.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030153.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030154.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030155.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030156.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030157.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030158.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030159.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030160.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030161.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030162.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030163.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030164.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030165.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030166.bam ; \
  /home/afarre/.local/necklace-necklace_v0.9/tools/bin/featureCounts \
  -T 32 --primary -p -t exon -g gene_id  --fraction -O -f \
  -a /home/afarre/Vem_necklace/mapped_reads/stringtie_merged_out.bam/stringtie_merged_bam_blocks.gtf \
  -o /home/afarre/Vem_necklace/mapped_reads/stringtie_merged_out.bam/block.counts \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030152.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030153.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030154.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030155.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030156.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030157.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030158.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030159.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030160.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030161.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030162.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030163.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030164.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030165.bam \
  /data/projects/punim0352/Vem_necklace/mapped_reads/DRR030166.bam 
  