#!/bin/bash
# STEP 9.1 -- allelic FIMO on the independent cis signals of >=2-signal genes
# (analysis #3). Per-SNP logic adapted from Crystal's proven FIMO script (bedtools
# slop +/-50bp, getfasta -name, explicit REF/ALT construction with allele checks,
# FIMO --thresh 1, keep motif hits overlapping the SNP, classify REF vs ALT change).
# Writes results/conditional_analysis/fimo/out/fimo_output_<snp>.txt per SNP.
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 06:00
#BSUB -q short
#BSUB -J fimo_signals
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step9.1_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step9.1_%J.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts
FIMODIR=${ROOT}/results/conditional_analysis/fimo
BEDS=${FIMODIR}/beds
OUT=${FIMODIR}/out
SIF=/share/pkg/containers/rstudio_example/r_4.5.2.sif

genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
genome_fa=/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme
# bedtools via module (as in step0_*); fimo via the local MEME install
module load bedtools 2>/dev/null || true
export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH
CONDA_PY=/home/shuo.shan-umw/miniconda3/envs/fastQTL/bin/python3   # pandas for the merge step
command -v bedtools >/dev/null || { echo "ERROR: bedtools not on PATH (adjust 'module load bedtools')"; exit 1; }
command -v fimo     >/dev/null || { echo "ERROR: fimo not on PATH (check MEME install path)"; exit 1; }

mkdir -p ${OUT}

# ---- 1) build per-SNP BEDs (R in Singularity) ------------------------------
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec ${SIF} Rscript ${SCRIPTDIR}/step9.1_build_signal_snp_beds.R ${ROOT}

# ---- 2) per-SNP FIMO loop --------------------------------------------------
while read -r snp; do
  [[ -z "$snp" ]] && continue
  snp_bed="${BEDS}/QTL_${snp}.bed"
  out="${OUT}/fimo_output_${snp}.txt"
  [[ ! -f "$snp_bed" ]] && { echo "no bed for $snp"; continue; }
  [[ -s "$out" ]] && { echo "skip (done) $snp"; continue; }   # idempotent

  wd=$(mktemp -d); pushd "$wd" >/dev/null
  cp "$snp_bed" SNP.bed

  bedtools slop -b 50 -i SNP.bed -g "$genome" > SNP_slop100.bed
  bedtools getfasta -fi "$genome_fa" -bed SNP_slop100.bed -fo SNP_slop100.fa -name

  # build REF and ALT fasta + variant intervals (stdlib python)
  python3 <<'PY'
with open("SNP.bed") as f: fields=f.readline().strip().split()
chrom,start,end,snp_id,ref,alt=fields[:6]; start=int(start); end=int(end); ref=ref.upper(); alt=alt.upper()
if "," in alt: raise ValueError("multi-allelic ALT not supported: "+alt)
with open("SNP_slop100.bed") as f: wchrom,wstart,wend=f.readline().strip().split()[:3]
wstart=int(wstart); wend=int(wend)
variant_idx0=start-wstart; variant_pos1=variant_idx0+1
lines=[x.strip() for x in open("SNP_slop100.fa") if x.strip()]; seq="".join(lines[1:]).upper()
if variant_idx0<0 or variant_idx0>=len(seq): raise ValueError("variant index out of range")
fasta_ref_span=seq[variant_idx0:variant_idx0+len(ref)]; fasta_alt_span=seq[variant_idx0:variant_idx0+len(alt)]
if fasta_ref_span==ref: anchor_len=len(ref)
elif fasta_alt_span==alt: anchor_len=len(alt)
else: raise ValueError(f"neither allele matches FASTA: bedREF={ref} span={fasta_ref_span} bedALT={alt} at {chrom}:{start}-{end}")
left=seq[:variant_idx0]; right=seq[variant_idx0+anchor_len:]
seq_ref=left+ref+right; seq_alt=left+alt+right
vtype="SNV" if len(ref)==1 and len(alt)==1 else ("deletion" if len(ref)>len(alt) else ("insertion" if len(ref)<len(alt) else "MNV"))
open("SNP_slop100_REF.fa","w").write(f">{snp_id}|REF|variant_pos={variant_pos1}|{vtype}\n{seq_ref}\n")
open("SNP_slop100_ALT.fa","w").write(f">{snp_id}|ALT|variant_pos={variant_pos1}|{vtype}\n{seq_alt}\n")
open("variant_interval_REF.txt","w").write(f"{variant_pos1}\t{variant_pos1+len(ref)-1}\n")
open("variant_interval_ALT.txt","w").write(f"{variant_pos1}\t{variant_pos1+len(alt)-1}\n")
PY

  for allele in REF ALT; do
    fimo --oc "results_${allele}" --thresh 1 --verbosity 1 "$memeF" "SNP_slop100_${allele}.fa"
    awk 'BEGIN{FS=OFS="\t"} !/^#/ && NR>1 && $1!=""' "results_${allele}/fimo.tsv" > "fimo_all_${allele}.tsv"
    vs=$(cut -f1 "variant_interval_${allele}.txt"); ve=$(cut -f2 "variant_interval_${allele}.txt")
    awk -v snp="$snp" -v allele="$allele" -v vs="$vs" -v ve="$ve" \
      'BEGIN{FS=OFS="\t"}{ if ($4<=ve && $5>=vs) print snp, allele, $0 }' \
      "fimo_all_${allele}.tsv" > "fimo_snp_overlap_${allele}.tsv"
  done

  # merge REF/ALT overlapping hits, classify per motif (pandas)
  "$CONDA_PY" <<'PY'
import pandas as pd, numpy as np
cols=["snp","allele","motif_id","motif_alt_id","sequence_name","start","stop","strand","score","pvalue","qvalue","matched_sequence"]
def rd(fn):
    try: return pd.read_csv(fn,sep="\t",header=None,names=cols)
    except pd.errors.EmptyDataError: return pd.DataFrame(columns=cols)
def best(df,s):
    if df.empty: return pd.DataFrame(columns=["motif_id",f"score_{s}",f"pvalue_{s}"])
    df=df.sort_values(["motif_id","pvalue"]).groupby("motif_id",as_index=False).first()
    return df.rename(columns={"score":f"score_{s}","pvalue":f"pvalue_{s}"})[["motif_id",f"score_{s}",f"pvalue_{s}"]]
ref=best(rd("fimo_snp_overlap_REF.tsv"),"REF"); alt=best(rd("fimo_snp_overlap_ALT.tsv"),"ALT")
m=ref.merge(alt,on="motif_id",how="outer")
m["present_REF"]=~m["score_REF"].isna(); m["present_ALT"]=~m["score_ALT"].isna()
def cls(r):
    if r["present_REF"] and not r["present_ALT"]: return "lost_in_ALT"
    if not r["present_REF"] and r["present_ALT"]: return "gained_in_ALT"
    if r["score_ALT"]>r["score_REF"]: return "stronger_in_ALT"
    if r["score_ALT"]<r["score_REF"]: return "weaker_in_ALT"
    return "unchanged"
m["neg_log10p_REF"]=-np.log10(m["pvalue_REF"].clip(lower=1e-300)); m["neg_log10p_ALT"]=-np.log10(m["pvalue_ALT"].clip(lower=1e-300))
m["delta_score_ALT_minus_REF"]=m["score_ALT"]-m["score_REF"]; m["delta_log10p_ALT_minus_REF"]=m["neg_log10p_ALT"]-m["neg_log10p_REF"]
m["priority_abs_delta_log10p"]=m["delta_log10p_ALT_minus_REF"].abs(); m["motif_change_class"]=m.apply(cls,axis=1)
m=m.sort_values(["motif_change_class","priority_abs_delta_log10p"],ascending=[True,False],na_position="last")
m.to_csv("fimo_REF_ALT_comparison.tsv",sep="\t",index=False)
PY

  if [[ -s fimo_REF_ALT_comparison.tsv ]]; then mv fimo_REF_ALT_comparison.tsv "$out"; echo "done $snp"; else echo "FAIL $snp"; fi
  popd >/dev/null; rm -rf "$wd"
done < "${FIMODIR}/signal_snp_list.txt"

echo "step9.1 FIMO done -> ${OUT}/fimo_output_<snp>.txt"; date
