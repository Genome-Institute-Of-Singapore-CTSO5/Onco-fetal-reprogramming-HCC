#1 - run on ionode, sort your bam files using samtools sort
/mnt/software/stow/samtools-1.7/bin/samtools sort -l 7 -m 2000M -t CB -O BAM -@ 16 -o /mnt/projects/seowjjw/singlecellDGLab/Regina/Fetal/XHL090/outs/cellsorted_possorted_genome_bam.bam /mnt/projects/seowjjw/singlecellDGLab/Regina/Fetal/XHL090/outs/possorted_genome_bam.bam

#2 - run velocyto, to generate .loom files 
/mnt/software/unstowable/miniconda3-4.6.14/envs/mamba/bin/velocyto run -b /mnt/projects/seowjjw/singlecellDGLab/Regina/Fetal/XHL089/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv -o /mnt/projects/seowjjw/singlecellDGLab/Regina/Fetal/XHL089/ -m /home/rmmwong/hg38_rmsk.gtf /mnt/projects/seowjjw/singlecellDGLab/Regina/Fetal/XHL089/outs/possorted_genome_bam.bam /mnt/projects/ycao/Data_2017/10x_scRNAseq/liver/refdata/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf

#3 - refer to jupyter notebook that merges loom files together

#4 - refer to jupyter notebook that plots using scVelo