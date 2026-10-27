#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=18000]
#BSUB -q long
#BSUB -W 12:00
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/tools/ABC-Enhancer-Gene-Prediction/logs/%J_download.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/tools/ABC-Enhancer-Gene-Prediction/logs/%J_download.err"

# ============================================================
#  Optional: fetch/convert 3D-contact maps for ABC.
#  You do NOT need any of this for the default (CONTACT_MODE=powerlaw).
#  Only run the section you want, then set CONTACT_MODE=hic + HIC_DIR
#  in the matching run_ABC_<CONTEXT>.sh scripts.
# ============================================================
set -euo pipefail

ABC_REPO='/pi/manuel.garber-umw/human/skin/eQTLs/tools/ABC-Enhancer-Gene-Prediction'
JUICER_JAR='/pi/manuel.garber-umw/human/skin/eQTLs/tools/juicer_tools.jar'   # download from Aiden lab if you don't have it
OUT_ROOT='/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/ABC/hic'
RES=5000
mkdir -p "${OUT_ROOT}"

# ------------------------------------------------------------
# make Gene bed
# ------------------------------------------------------------
GTF=/pi/manuel.garber-umw/home/yuqing.wang-umw/inDrop/Skin/diablo/build/hg38_v34/gencode.v34.primary_assembly.annotation.gtf

awk -F'\t' '$3=="gene" && $1 ~ /^chr([0-9]+|[XY])$/ {
    name=$9; sub(/.*gene_name "/,"",name); sub(/".*/,"",name);
    print $1"\t"($4-1)"\t"$5"\t"name"\t0\t"$7
}' "$GTF" | sort -k1,1 -k2,2n \
  > /pi/manuel.garber-umw/human/skin/eQTLs/annotation/gencode_v34_gene_bounds.bed

# ------------------------------------------------------------
# (A) FIBROBLAST — IMR90 in situ Hi-C (Rao 2014, GSE63525). Genome-wide, drop-in.
#     Real GEO supplementary .hic; then dump to ABC's per-chrom KR-normalized dir.
# ------------------------------------------------------------
imr90() {
  local D="${OUT_ROOT}/IMR90"; mkdir -p "$D"
  wget -c -O "${D}/GSE63525_IMR90_combined_30.hic" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_IMR90_combined_30.hic"
  python "${ABC_REPO}/src/juicebox_dump.py" \
    --hic_file "${D}/GSE63525_IMR90_combined_30.hic" \
    --juicebox "java -jar ${JUICER_JAR}" \
    --outdir "${OUT_ROOT}/FRB" --resolution ${RES}
  echo "IMR90 -> ${OUT_ROOT}/FRB  (set HIC_DIR to this in run_ABC_FRB_*.sh; note: fetal lung, a PROXY for dermal FRB)"
}

# ------------------------------------------------------------
# (B) AVERAGE Hi-C (Nasser 2021) — cell-type-agnostic, near matched-Hi-C performance.
#     Distributed by the ABC/Engreitz team; see ABC ReadTheDocs "average Hi-C".
#     Put the averaged bedpe/dir here and point every HIC_DIR at it.
# ------------------------------------------------------------
average_hic() {
  echo "Fetch the ABC average Hi-C per the docs: https://abc-enhancer-gene-prediction.readthedocs.io"
  echo "Then set CONTACT_MODE=hic and HIC_DIR=<avg dir> in all run_ABC_*.sh (use for KRT/MEL which lack matched Hi-C)."
}

# ------------------------------------------------------------
# (C) KERATINOCYTE — H3K27ac HiChIP (best cell-type match for KRT).
#     TODO: get GEO/accession from the paper data-availability, then either
#       (i) call loops with FitHiChIP and convert to an ABC HiCdir, or
#      (ii) use the ABC HiChIP-as-contact path.
#     Refs: keratinocyte capture Hi-C JACI S009167492031410X ; HaCaT psoriasis PMID 32366252.
# ------------------------------------------------------------
keratinocyte_hichip() { echo "See README section 'Keratinocyte HiChIP' — accession + FitHiChIP steps."; }

# choose what to run:
#   bash download_contacts.sh imr90
#   bash download_contacts.sh average_hic
"${1:-}" || { echo "usage: bash download_contacts.sh [imr90|average_hic|keratinocyte_hichip]"; exit 0; }
