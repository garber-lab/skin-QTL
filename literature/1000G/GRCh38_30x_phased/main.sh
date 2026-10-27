#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 8:00
#BSUB -o /pi/manuel.garber-umw/human/skin/eQTLs/literature/1000G/GRCh38_30x_phased/logs/main.%I.out
#BSUB -e /pi/manuel.garber-umw/human/skin/eQTLs/literature/1000G/GRCh38_30x_phased/logs/main.%I.err

# scripts to download the 1000G phased VCF file and metadata
ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/literature/1000G/GRCh38_30x_phased
mkdir -p $ROOT/vcf $ROOT/metadata $ROOT/logs

# download phased vcf, one chromosome per job
cd ${ROOT}
bsub < download_1000g_hg38_phased.lsf.sh


# download metadata
mkdir -p "$ROOT/metadata"

wget -c -P "$ROOT/metadata" \
  "$BASE/1kGP.3202_samples.pedigree_info.txt"

wget -c -P "$ROOT/metadata" \
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"


# validate after download
#VCF=$ROOT/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz
#
#bcftools query -l "$VCF" | wc -l
#tabix -l "$VCF" | head
#bcftools view -H "$VCF" chr1:1000000-1010000 | head
