#!/bin/bash
export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages

base_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5
data_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/data
shap_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/variant_shap/rs5765171
shap_variant_input=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/variant_shap/rs5765171/shap_variant_input.txt
variant_scorer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src
HOCOMOCO=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme
VARIANT_ID=rs5765171
GENE=PHF21B
CONTEXT=MEL_IFNG
FLANK_LOGO=60

echo "================================================"; date
echo "Variant SHAP: ${GENE} ${VARIANT_ID} in ${CONTEXT}"
echo "================================================"

# ------------------------------------------------------------
# STEP 2: per-fold SHAP, counts and profile
# ------------------------------------------------------------
for fold_n in 0 1 2 3 4; do
    model=${base_dir}/fold${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5
    fold_shap_dir=${shap_dir}/fold${fold_n}
    mkdir -p ${fold_shap_dir}

    if [ ! -f "${model}" ]; then
        echo "WARNING: fold ${fold_n} model not found, skipping..."; continue
    fi

    for shap_type in counts profile; do
        if [ -f "${fold_shap_dir}/shap.variant_shap.${shap_type}.h5" ]; then
            echo "fold ${fold_n} ${shap_type} SHAP exists, skipping..."
            continue
        fi
        echo "--- fold ${fold_n}: ${shap_type} SHAP ---"; date
        python ${variant_scorer}/variant_shap.py \
            --list        ${shap_variant_input} \
            --genome      ${data_dir}/downloads/hg38.fa \
            --model       ${model} \
            --out_prefix  ${fold_shap_dir}/shap \
            --chrom_sizes ${data_dir}/downloads/hg38.chrom.sizes \
            --schema      original \
            --shap_type   ${shap_type}
        echo "--- fold ${fold_n} ${shap_type} done ---"; date
    done
done

# ------------------------------------------------------------
# STEP 3: average across folds
#   h5 layout from variant_shap.py:
#     raw/seq            (2, 4, 2114)  one-hot, identical across folds
#     shap/seq           (2, 4, 2114)  hypothetical importance
#     projected_shap/seq (2, 4, 2114)  hyp * one-hot
#     alleles            (2,)          0 = allele1 (REF), 1 = allele2 (ALT)
# ------------------------------------------------------------
merged_dir=${shap_dir}/merged
mkdir -p ${merged_dir}

for shap_type in counts profile; do
    if [ -f "${merged_dir}/shap_${shap_type}.npz" ] && [ -f "${merged_dir}/ohe_${shap_type}.npz" ]; then
        echo "merged ${shap_type} npz exists, skipping..."; continue
    fi
    echo "--- averaging ${shap_type} across folds ---"; date

python3 << PYEOF
import h5py, numpy as np, sys
try:
    import hdf5plugin
except ImportError:
    pass

shap_dir   = "${shap_dir}"
shap_type  = "${shap_type}"
merged_dir = "${merged_dir}"

ohe_ref = shap_acc = proj_acc = None
n_folds = 0
for fold_n in range(5):
    p = f"{shap_dir}/fold{fold_n}/shap.variant_shap.{shap_type}.h5"
    try:
        with h5py.File(p, "r") as f:
            ohe  = f["raw/seq"][:]
            shap = f["shap/seq"][:]
            proj = f["projected_shap/seq"][:]
    except Exception as e:
        print(f"  fold{fold_n}: skipping ({e})"); continue
    if ohe_ref is None:
        ohe_ref, shap_acc, proj_acc = ohe, shap.copy(), proj.copy()
    else:
        shap_acc += shap; proj_acc += proj
    n_folds += 1
    print(f"  fold{fold_n}: shap range [{shap.min():.4f}, {shap.max():.4f}]")

if n_folds == 0:
    print(f"ERROR: no {shap_type} h5 found"); sys.exit(1)

np.savez_compressed(f"{merged_dir}/ohe_{shap_type}.npz",  ohe_ref)
np.savez_compressed(f"{merged_dir}/shap_{shap_type}.npz", shap_acc / n_folds)
np.savez_compressed(f"{merged_dir}/proj_{shap_type}.npz", proj_acc / n_folds)
print(f"averaged over {n_folds} folds -> {merged_dir}")
PYEOF

done

# ------------------------------------------------------------
# STEP 4: modisco -- seqlet extraction + HOCOMOCO matching.
# See the header: with n=2 sequences this is motif ANNOTATION,
# not de novo discovery.
# ------------------------------------------------------------
for shap_type in counts profile; do
    modisco_out=${shap_dir}/merged/modisco_${shap_type}
    if [ -f "${modisco_out}/modisco_results.h5" ]; then
        echo "${shap_type} modisco exists, skipping..."; continue
    fi
    ohe_npz=${merged_dir}/ohe_${shap_type}.npz
    shap_npz=${merged_dir}/shap_${shap_type}.npz
    [ -f "${ohe_npz}" ] && [ -f "${shap_npz}" ] || { echo "missing npz, skipping"; continue; }

    mkdir -p ${modisco_out}/report
    echo "--- modisco motifs: ${shap_type} ---"; date
    env -u PYTHONPATH modisco motifs \
        -s ${ohe_npz} \
        -a ${shap_npz} \
        -n 2 \
        -w 200 \
        -o ${modisco_out}/modisco_results.h5

    echo "--- modisco report: ${shap_type} ---"; date
    env -u PYTHONPATH modisco report \
        -i ${modisco_out}/modisco_results.h5 \
        -o ${modisco_out}/report/ \
        -s ${modisco_out}/report/ \
        -m ${HOCOMOCO}
done

# ------------------------------------------------------------
# STEP 5: sequence logos, projected SHAP, REF vs ALT
# variant sits at index 1057 of the 2114 bp input
# ------------------------------------------------------------
for shap_type in counts profile; do
    out_pdf=${merged_dir}/logos_${shap_type}.pdf
    if [ -f "${out_pdf}" ]; then echo "logos_${shap_type} exists, skipping..."; continue; fi

env -u PYTHONPATH python3 << PYEOF
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker, pandas as pd

merged_dir = "${merged_dir}"
shap_type  = "${shap_type}"
flank      = int("${FLANK_LOGO}")
gene       = "${GENE}"
vid        = "${VARIANT_ID}"
ctx        = "${CONTEXT}"

proj = np.load(f"{merged_dir}/proj_{shap_type}.npz")["arr_0"]   # (2, 4, 2114)
L      = proj.shape[2]
center = L // 2
lo, hi = center - flank, center + flank + 1

fig, axes = plt.subplots(2, 1, figsize=(16, 5), sharex=True, sharey=True)
for i, lab in enumerate(["allele1 (REF)", "allele2 (ALT)"]):
    df = pd.DataFrame(proj[i, :, lo:hi].T, columns=["A", "C", "G", "T"])
    df.index = np.arange(lo - center, hi - center)
    logomaker.Logo(df, ax=axes[i])
    axes[i].axvline(0, ls="--", lw=1, color="red")
    axes[i].set_ylabel(lab)
axes[0].set_title(f"{gene}  {vid}  {ctx}  --  projected SHAP ({shap_type})")
axes[1].set_xlabel(f"position relative to variant (bp), index {center} of {L}")
plt.tight_layout()
plt.savefig(f"{merged_dir}/logos_{shap_type}.pdf")
plt.savefig(f"{merged_dir}/logos_{shap_type}.png", dpi=150)
print(f"wrote logos_{shap_type}.pdf/.png")
PYEOF

done

echo "================================================"; date
echo "DONE  ${GENE} ${VARIANT_ID} ${CONTEXT}"
echo "  per-fold h5   : ${shap_dir}/fold{0..4}/"
echo "  merged npz    : ${shap_dir}/merged/"
echo "  modisco report: ${shap_dir}/merged/modisco_{counts,profile}/report/"
echo "  logos         : ${shap_dir}/merged/logos_{counts,profile}.pdf"
echo "================================================"
