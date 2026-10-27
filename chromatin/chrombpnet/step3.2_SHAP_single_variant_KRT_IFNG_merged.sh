#!/bin/bash
#BSUB -J chrombpnet_shap_KRT_IFNG
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -W 2:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step3.2_shap_KRT_IFNG.%J.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step3.2_shap_KRT_IFNG.%J.out"

# ============================================================
# VARIANT SHAP — single targeted variant
#
# Purpose: compute per-position DeepLIFT/SHAP importance scores
# for both alleles of a specific variant, to identify which TF
# motif(s) are disrupted.
#
# Target: rs11204676 (CTSS eQTL, KRT_IFNG, top ChromBPNet hit)
# Variant ID in scored list:
#   rs11204676_chr16_1930267_G_A   (G = REF/allele1, A = ALT/allele2)
#
# All outputs are written under:
#   ${base_shap_dir}/${VARIANT_ID}/
# ============================================================

# ============================================================
# TARGET VARIANT
# Edit VARIANT_ID to target a different SNP. Use rsID prefix.
# Must match the start of a variant_id in variant_list_original_schema.txt
# ============================================================
VARIANT_ID="rs11177432"

# ============================================================
# PATHS — inherited from step3 conventions
# ============================================================
BIAS_FACTOR=0.5
sample=KRT_IFN
prefix=ATAC_${sample}_merged_alldonors

base_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}/biasFactor${BIAS_FACTOR}"
data_dir="${base_dir}/data"
score_dir="${base_dir}/variant_scores"
base_shap_dir="${base_dir}/variant_shap"
shap_dir="${base_shap_dir}/${VARIANT_ID}"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

variant_scorer="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src"
conda_site_pkgs="/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages"

HOCOMOCO="/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme"

# ============================================================
# STEP 1: EXTRACT SINGLE-VARIANT INPUT FILE
# Reuses the already-prepared and allele-validated input from step3.
# ============================================================
echo "--- Extracting variant ${VARIANT_ID} ---"; date

mkdir -p "${shap_dir}"

shap_variant_input="${shap_dir}/shap_variant_input.txt"

grep ${VARIANT_ID}_ "${score_dir}/variant_list_original_schema.txt" > "${shap_variant_input}"

n_found=$(wc -l < "${shap_variant_input}")
echo "Variants found: ${n_found}"

if [ "${n_found}" -eq 0 ]; then
    echo "ERROR: ${VARIANT_ID} not found in ${score_dir}/variant_list_original_schema.txt"
    echo "Check the variant ID or re-run step3 variant prep."
    exit 1
fi

cat "${shap_variant_input}"

# ============================================================
# STEP 2: RUN SHAP INSIDE SINGULARITY, ONCE PER FOLD
# No --bias arg (variant_shap.py operates on the nobias model only)
# No --num_shuf (not applicable to SHAP)
# No --peaks   (not a valid argument for variant_shap.py)
# ============================================================
rm -f "${base_dir}/run_inside_singularity_shap.sh"
cat > "${base_dir}/run_inside_singularity_shap.sh" << EOF
#!/bin/bash
export PYTHONPATH=${conda_site_pkgs}

base_dir=${base_dir}
data_dir=${data_dir}
shap_dir=${shap_dir}
prefix=${prefix}
shap_variant_input=${shap_variant_input}
variant_scorer=${variant_scorer}
HOCOMOCO=${HOCOMOCO}

echo "================================================"; date
echo "Variant SHAP for ${VARIANT_ID}"
echo "================================================"

# ============================================================
# STEP 2: RUN SHAP PER FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=\${base_dir}/fold\${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5
    fold_shap_dir=\${shap_dir}/fold\${fold_n}
    mkdir -p \${fold_shap_dir}

    if [ ! -f "\${model}" ]; then
        echo "WARNING: fold \${fold_n} model not found, skipping..."; continue
    fi

    if [ -f "\${fold_shap_dir}/shap.variant_shap.counts.h5" ]; then
        echo "fold \${fold_n} counts SHAP already exists, skipping..."
    else
        echo "--- fold \${fold_n}: counts SHAP ---"; date
        python \${variant_scorer}/variant_shap.py \\
            --list       \${shap_variant_input} \\
            --genome     \${data_dir}/downloads/hg38.fa \\
            --model      \${model} \\
            --out_prefix \${fold_shap_dir}/shap \\
            --chrom_sizes \${data_dir}/downloads/hg38.chrom.sizes \\
            --schema     original \\
            --shap_type  counts
        echo "--- fold \${fold_n} counts SHAP done ---"; date
    fi

    if [ -f "\${fold_shap_dir}/shap.variant_shap.profile.h5" ]; then
        echo "fold \${fold_n} profile SHAP already exists, skipping..."
    else
        echo "--- fold \${fold_n}: profile SHAP ---"; date
        python \${variant_scorer}/variant_shap.py \\
            --list       \${shap_variant_input} \\
            --genome     \${data_dir}/downloads/hg38.fa \\
            --model      \${model} \\
            --out_prefix \${fold_shap_dir}/shap \\
            --chrom_sizes \${data_dir}/downloads/hg38.chrom.sizes \\
            --schema     original \\
            --shap_type  profile
        echo "--- fold \${fold_n} profile SHAP done ---"; date
    fi

done

echo "================================================"
echo "SHAP outputs written to: \${shap_dir}/fold{0..4}/"
ls -lh \${shap_dir}/fold0/ 2>/dev/null || echo "  (fold0 output not found)"
echo "================================================"

# ============================================================
# STEP 3: AVERAGE SHAP TRACKS ACROSS FOLDS
#
# h5 structure from variant_shap.py:
#   raw/seq            (2, 4, 2114)  one-hot — sequence-derived,
#                                    identical across folds, take from fold0
#   shap/seq           (2, 4, 2114)  hypothetical importance scores
#   projected_shap/seq (2, 4, 2114)  hyp * one-hot (actual importances)
#   alleles            (2,)          0=allele1, 1=allele2
#
# Average shap/seq and projected_shap/seq across available folds.
# Write merged npz files ready for modisco and logo plotting.
# ============================================================
echo "================================================"; date
echo "Averaging SHAP tracks across folds"
echo "================================================"

merged_dir=\${shap_dir}/merged
mkdir -p \${merged_dir}

for shap_type in counts profile; do

    if [ -f "\${merged_dir}/ohe_\${shap_type}.npz" ] && [ -f "\${merged_dir}/shap_\${shap_type}.npz" ]; then
        echo "merged \${shap_type} npz already exists, skipping..."
        continue
    fi

    echo "--- averaging \${shap_type} across folds ---"; date

    python3 << PYEOF
import h5py, numpy as np, sys
try:
    import hdf5plugin  # registers Blosc/LZ4 codecs used by variant_shap.py
except ImportError:
    pass

shap_dir   = "\${shap_dir}"
shap_type  = "\${shap_type}"
merged_dir = "\${merged_dir}"

ohe_ref  = None
shap_acc = None
proj_acc = None
n_folds  = 0

for fold_n in range(5):
    h5_path = f"{shap_dir}/fold{fold_n}/shap.variant_shap.{shap_type}.h5"
    try:
        with h5py.File(h5_path, "r") as f:
            ohe  = f["raw/seq"][:]             # (2, 4, 2114)
            shap = f["shap/seq"][:]            # (2, 4, 2114)
            proj = f["projected_shap/seq"][:]  # (2, 4, 2114)
    except Exception as e:
        print(f"  fold{fold_n}: skipping ({e})")
        continue

    if ohe_ref is None:
        ohe_ref  = ohe
        shap_acc = shap.copy()
        proj_acc = proj.copy()
    else:
        shap_acc += shap
        proj_acc += proj
    n_folds += 1
    print(f"  fold{fold_n}: loaded  shap range [{shap.min():.4f}, {shap.max():.4f}]")

if n_folds == 0:
    print(f"ERROR: no {shap_type} h5 files found — did SHAP complete?")
    sys.exit(1)

shap_mean = shap_acc / n_folds
proj_mean = proj_acc / n_folds

np.savez_compressed(f"{merged_dir}/ohe_{shap_type}.npz",  ohe_ref)
np.savez_compressed(f"{merged_dir}/shap_{shap_type}.npz", shap_mean)
np.savez_compressed(f"{merged_dir}/proj_{shap_type}.npz", proj_mean)

print(f"Averaged over {n_folds} folds")
print(f"  ohe:  {ohe_ref.shape}")
print(f"  shap: {shap_mean.shape}  range [{shap_mean.min():.4f}, {shap_mean.max():.4f}]")
print(f"  proj: {proj_mean.shape}  range [{proj_mean.min():.4f}, {proj_mean.max():.4f}]")
print(f"  Written to {merged_dir}/")
PYEOF

    echo "--- \${shap_type} averaging done ---"; date
done

# ============================================================
# STEP 4: RUN MODISCO-LITE ON AVERAGED TRACKS
#
# Both alleles (n=2) are passed together so modisco can extract
# one seqlet per allele and the report shows both side-by-side
# against the HOCOMOCOv12 motif database.
#
# -n 2  : max seqlets — one per allele
# -w 200: +/-100bp window around each seqlet center
#
# env -u PYTHONPATH: prevents conda site-packages from bleeding
# into the container's modisco installation
# ============================================================
echo "================================================"; date
echo "Running modisco-lite on averaged SHAP tracks"
echo "================================================"

for shap_type in counts profile; do
    modisco_out=\${shap_dir}/merged/modisco_\${shap_type}

    if [ -f "\${modisco_out}/modisco_results.h5" ]; then
        echo "\${shap_type} modisco already exists, skipping..."; continue
    fi

    ohe_npz=\${merged_dir}/ohe_\${shap_type}.npz
    shap_npz=\${merged_dir}/shap_\${shap_type}.npz

    if [ ! -f "\${ohe_npz}" ] || [ ! -f "\${shap_npz}" ]; then
        echo "ERROR: merged npz for \${shap_type} not found, skipping modisco"
        continue
    fi

    mkdir -p \${modisco_out}/report

    echo "--- modisco motifs: \${shap_type} ---"; date
    env -u PYTHONPATH modisco motifs \\
        -s \${ohe_npz} \\
        -a \${shap_npz} \\
        -n 2 \\
        -w 200 \\
        -o \${modisco_out}/modisco_results.h5

    echo "--- modisco report: \${shap_type} ---"; date
    env -u PYTHONPATH modisco report \\
        -i \${modisco_out}/modisco_results.h5 \\
        -o \${modisco_out}/report/ \\
        -s \${modisco_out}/report/ \\
        -m \${HOCOMOCO}

    echo "--- modisco \${shap_type} done ---"; date
done

# ============================================================
# STEP 5: SEQUENCE LOGO PLOTS — projected SHAP per allele
#
# For each shap_type (counts, profile), plot a +/-60bp window
# around the variant center (index 1057 of 2114) showing the
# projected importance scores (hyp * one-hot) for both alleles
# side by side.  Row 0 = allele1 (REF), row 1 = allele2 (ALT).
# The variant position is marked with a dashed vertical line.
#
# env -u PYTHONPATH: prevents conda pandas/matplotlib from
# bleeding into the container's logomaker installation
#
# Output: merged/logos_{counts,profile}.pdf + .png
# ============================================================
echo "================================================"; date
echo "Plotting sequence logos"
echo "================================================"

env -u PYTHONPATH python3 << PYEOF2
import numpy as np, os, sys
try:
    import hdf5plugin
except ImportError:
    pass

shap_dir   = "\${shap_dir}"
merged_dir = "\${merged_dir}"
variant_id = "${VARIANT_ID}"

center   = 1057   # center of 2114bp window (0-indexed)
half_win = 60     # bp on each side to plot
bases    = list("ACGT")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import logomaker
import pandas as pd

def make_logo_df(proj_4xL, pos_start, pos_end):
    window = proj_4xL[:, pos_start:pos_end].T   # (L, 4)
    df = pd.DataFrame(window, columns=bases)
    df.index = range(pos_start - center, pos_end - center)
    return df

for shap_type in ["counts", "profile"]:
    proj_path = os.path.join(merged_dir, f"proj_{shap_type}.npz")
    if not os.path.exists(proj_path):
        print(f"  {shap_type}: proj npz not found, skipping plot")
        continue

    proj     = np.load(proj_path)["arr_0"]   # (2, 4, 2114)
    proj_ref = proj[0]                        # allele1 REF
    proj_alt = proj[1]                        # allele2 ALT

    pos_start = center - half_win
    pos_end   = center + half_win + 1

    df_ref = make_logo_df(proj_ref, pos_start, pos_end)
    df_alt = make_logo_df(proj_alt, pos_start, pos_end)

    # shared y-axis limits across both panels
    ymax = max(df_ref.clip(lower=0).sum(axis=1).max(),
               df_alt.clip(lower=0).sum(axis=1).max()) * 1.15
    ymin = min(df_ref.clip(upper=0).sum(axis=1).min(),
               df_alt.clip(upper=0).sum(axis=1).min()) * 1.15
    if ymax == 0: ymax =  0.01
    if ymin == 0: ymin = -0.01

    fig, axes = plt.subplots(2, 1, figsize=(14, 5), sharex=True,
                             gridspec_kw={"hspace": 0.08})

    allele_labels = ["allele1  REF  (G)", "allele2  ALT  (A)"]
    for ax, df, label in zip(axes, [df_ref, df_alt], allele_labels):
        logomaker.Logo(df, ax=ax, color_scheme="classic", baseline_width=0.5)
        ax.axvline(0, color="black", linewidth=1.2, linestyle="--", alpha=0.7)
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel("importance\nscore", fontsize=9)
        ax.text(0.01, 0.88, label, transform=ax.transAxes, fontsize=10, va="top")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

    axes[1].set_xlabel("position relative to variant (bp)", fontsize=10)
    fig.suptitle(
        f"{variant_id}  |  ChromBPNet projected SHAP  |  {shap_type}",
        fontsize=11, y=1.01
    )

    for ext in ["pdf", "png"]:
        out = os.path.join(merged_dir, f"logos_{shap_type}.{ext}")
        fig.savefig(out, bbox_inches="tight", dpi=150)
        print(f"  saved: {out}")
    plt.close(fig)

print("Logo plots done.")
PYEOF2

echo "--- logo plots done ---"; date

echo "================================================"
echo "All done. Summary of outputs:"
echo "  Per-fold SHAP h5s  : \${shap_dir}/fold{0..4}/"
echo "  Averaged npz       : \${shap_dir}/merged/"
echo "  modisco results    : \${shap_dir}/merged/modisco_{counts,profile}/"
echo "  HTML reports       : \${shap_dir}/merged/modisco_{counts,profile}/report/"
echo "  Sequence logos     : \${shap_dir}/merged/logos_{counts,profile}.{pdf,png}"
echo "================================================"
date
EOF

# ============================================================
# LAUNCH INSIDE SINGULARITY
# ============================================================
echo "Waiting for filesystem sync..."; date; sleep 5

module load conda
source activate fastQTL

singularity exec --nv \
    --env PYTHONPATH=${conda_site_pkgs} \
    ${sif} bash ${base_dir}/run_inside_singularity_shap.sh

conda deactivate
