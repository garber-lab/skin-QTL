#!/bin/bash
# eQTL integration: align eQTL genotypes (old VCF) to the phased hg38 VCF, one chromosome per array task.
#
#   bsub -J "align[1-22]" ... < align_alleles.sh        (see the submit command in the plan / chat)
#
# Writes $OUT/chr<N>.align.tsv.gz (see align_alleles.py for the columns) and $OUT/chr<N>.log.
# After all 22 tasks: zcat chr*.align.tsv.gz | awk 'NR==1 || $1!="chrom"' | gzip > all.align.tsv.gz
set -euo pipefail
module load bcftools

E=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome
SIF=$E/../chrombpnet/chrombpnet_latest.sif
G=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output
OLD=$G/skineQTL_imputed_hg38_filtered.vcf.gz                                 # used by the eQTL pipeline
CHR=chr${LSB_JOBINDEX}
NEW=$G/phased/${CHR}_imputed_combined_hg38_rsID_filtered_phased.vcf.gz       # used by ChromBPNet
OUT=$E/eqtl_integration/align
TMP=$OUT/tmp_${CHR}
mkdir -p $TMP

# 1) genotypes of this chromosome from both VCFs (ID column = rsID in both)
FMT='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'
bcftools query -l $OLD > $TMP/old_samples.txt
bcftools query -l $NEW > $TMP/new_samples.txt
bcftools query -r $CHR -f "$FMT" $OLD | gzip > $TMP/old.tsv.gz
bcftools query -r $CHR -f "$FMT" $NEW | gzip > $TMP/new.tsv.gz

# 2) match SNPs, compare alleles and genotypes, decide the beta orientation
singularity exec $SIF python $E/eqtl_integration/scripts/align_alleles.py \
    $TMP/old.tsv.gz $TMP/old_samples.txt $TMP/new.tsv.gz $TMP/new_samples.txt \
    $OUT/${CHR}.align.tsv.gz > $OUT/${CHR}.log 2>&1

rm -r $TMP
