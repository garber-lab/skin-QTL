#!/usr/bin/env bash
#BSUB -J "dl_1000g[1-22]"
#BSUB -n 1
#BSUB -R "rusage[mem=2000]"
#BSUB -W 08:00
#BSUB -o /pi/manuel.garber-umw/human/skin/eQTLs/literature/1000G/GRCh38_30x_phased/logs/dl_1000g.%I.out
#BSUB -e /pi/manuel.garber-umw/human/skin/eQTLs/literature/1000G/GRCh38_30x_phased/logs/dl_1000g.%I.err

set -euo pipefail

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/literature/1000G/GRCh38_30x_phased
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased

CHR=chr${LSB_JOBINDEX}
FILE=CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz

mkdir -p "$ROOT/vcf" "$ROOT/logs"

wget -c -P "$ROOT/vcf" "$BASE/$FILE"
wget -c -P "$ROOT/vcf" "$BASE/$FILE.tbi"
