#!/bin/bash
#BSUB -J chrombpnet-b0.7
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5000]"
#BSUB -q gpu
#BSUB -W 720:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS_dedupped_bias_threshold_factor_0.7/log/pipeline.%J.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS_dedupped_bias_threshold_factor_0.7/log/pipeline.%J.out"

# preprocess ChromBPNet tutorial data
dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS_dedupped_bias_threshold_factor_0.7"
cd ${dir}
mkdir -p ${dir}/data
mkdir -p ${dir}/data/downloads

blacklist=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/blacklist.bed.gz
chromsize=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes
cp ${blacklist} ${dir}/data/downloads
cp ${chromsize} ${dir}/data/downloads
cp /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa ${dir}/data/downloads
cp /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa.fai ${dir}/data/downloads

### getting bam files for sample
module load samtools
prefix=ATAC_F111_KRT_PBS_S1
origbamF=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped/${prefix}.nodup.bam
origbamIdxF=${origbamF}.bai
echo "sorting bam now";date
cp ${origbamF} ${dir}/data/downloads/merged.bam
#samtools sort -@4 ${origbamF} -o ${dir}/data/downloads/merged.bam
echo "sorting done, indexing now"; date
samtools index ${dir}/data/downloads/merged.bam
echo "indexing done, processing peak data";date

### getting peak data (they recommend using relaxed peak-calls)
# MACS2 with p-val = 0.01
peakF=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2/${prefix}.clipped.narrowPeak.gz
zcat ${peakF} | sort -k1,1 -k2,2n | gzip > ${dir}/data/downloads/orig_peak.bed.gz 
# ensure that the peaks don't overlap blacklist
module load bedtools
bedtools slop -i ${dir}/data/downloads/blacklist.bed.gz -g ${dir}/data/downloads/hg38.chrom.sizes -b 1057 > ${dir}/data/downloads/temp.bed
bedtools intersect -v -a ${dir}/data/downloads/orig_peak.bed.gz -b ${dir}/data/downloads/temp.bed > ${dir}/data/downloads/peaks_no_blacklist.bed
echo "done with preprocessing step";date

### start singularity image to run chrombpnet
#### activate singularity
singularity exec --nv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif \
        bash -c '
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS_dedupped_bias_threshold_factor_0.7
### define train, validation and test chr splits
head -n 24 ${dir}/data/downloads/hg38.chrom.sizes > ${dir}/data/downloads/hg38.chrom.subset.sizes
mkdir -p ${dir}/data/splits
chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ${dir}/data/splits/fold_0
#chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr2 chr8 chr9 chr16 -vcr chr12 chr17 -op ${dir}/data/splits/fold_1
#chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr4 chr11 chr12 chr15 chrY -vcr chr22 chr7 -op ${dir}/data/splits/fold_2
#chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr5 chr10 chr14 chr18 chr20 chr22 -vcr chr6 chr21 -op ${dir}/data/splits/fold_3
#chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr7 chr13 chr17 chr19 chr21 chrX -vcr chr10 chr18 -op ${dir}/data/splits/fold_4 

### generate non-peaks (background regions)
hg38fasta=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa
cp /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa ${dir}/data/downloads/hg38.fa
cp /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa.fai ${dir}/data/downloads/hg38.fa.fai
chrombpnet prep nonpeaks -g ${hg38fasta} -p ${dir}/data/downloads/peaks_no_blacklist.bed -c ${dir}/data/downloads/hg38.chrom.sizes -fl ${dir}/data/splits/fold_0.json -br ${dir}/data/downloads/blacklist.bed.gz -o ${dir}/data/output
echo "done";date

### train bias model
mkdir -p ${dir}/bias_model
mkdir -p ${dir}/data/output
prefix=ATAC_F111_KRT_PBS_S1

### ------ train a bias model
echo "training bias model now"; date
chrombpnet bias pipeline \
        -ibam ${dir}/data/downloads/merged.bam \
        -d "ATAC" \
        -g ${dir}/data/downloads/hg38.fa \
        -c ${dir}/data/downloads/hg38.chrom.sizes \
        -p ${dir}/data/downloads/peaks_no_blacklist.bed \
        -n ${dir}/data/output_negatives.bed \
        -fl ${dir}/data/splits/fold_0.json \
        -b 0.6 \
        -o ${dir}/bias_model/ \
        -s 42 \
        -fp ${prefix}_fold0
echo "done";date


### ------ train bias-factorized ChromBPNet model
echo "training bias-factorized model using a pre-trained bias model now"; date
chrombpnet pipeline \
        -ibam ${dir}/data/downloads/merged.bam \
        -d "ATAC" \
        -g ${dir}/data/downloads/hg38.fa \
        -c ${dir}/data/downloads/hg38.chrom.sizes \
        -p ${dir}/data/downloads/peaks_no_blacklist.bed \
        -n ${dir}/data/output_negatives.bed \
        -fl ${dir}/data/splits/fold_0.json \
        -b ${dir}/bias_model/models/ATAC_F111_KRT_PBS_S1_fold0_bias.h5 \
        -o ${dir}/chrombpnet_model/
'

##---------- archive: interactive session
## start the gpu session
#bsub -Is -q gpu -W 8:00 -R rusage[mem=100G] /bin/bash
## start the singularity image
#singularity shell --nv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif
#
#### define train, validation and test chr splits
#dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS
#head -n 24 ${dir}/data/downloads/hg38.chrom.sizes > ${dir}/data/downloads/hg38.chrom.subset.sizes
#mkdir -p ${dir}/data/splits
#chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ${dir}/data/splits/fold_0
##chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr2 chr8 chr9 chr16 -vcr chr12 chr17 -op ${dir}/data/splits/fold_1
##chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr4 chr11 chr12 chr15 chrY -vcr chr22 chr7 -op ${dir}/data/splits/fold_2
##chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr5 chr10 chr14 chr18 chr20 chr22 -vcr chr6 chr21 -op ${dir}/data/splits/fold_3
##chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr7 chr13 chr17 chr19 chr21 chrX -vcr chr10 chr18 -op ${dir}/data/splits/fold_4 
#
#### generate non-peaks (background regions)
#hg38fasta=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa
#chrombpnet prep nonpeaks -g ${hg38fasta} -p ${dir}/data/downloads/peaks_no_blacklist.bed -c ${dir}/data/downloads/hg38.chrom.sizes -fl ${dir}/data/splits/fold_0.json -br ${dir}/data/downloads/blacklist.bed.gz -o ${dir}/data/output
#
