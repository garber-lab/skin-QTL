#!/bin/bash
#BSUB -J chrombpnet_footprints_KRT_PBS
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -q gpu
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step4_footprints_KRT_PBS_merged.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step4_footprints_KRT_PBS_merged.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE — only touch these when changing runs
# ============================================================
BIAS_FACTOR=0.5
sample=KRT_PBS
prefix=ATAC_${sample}_merged_alldonors

base_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}/biasFactor${BIAS_FACTOR}"
data_dir="${base_dir}/data"
footprint_dir="${base_dir}/footprints"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

# TF list — one TF name per line
tf_list_input="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/tf_list.txt"

# HOCOMOCO motif database
meme_file="/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme"

mkdir -p ${footprint_dir}

# ============================================================
# COPY TF LIST TO WORKING DIRECTORY FOR PROVENANCE
# ============================================================
cp ${tf_list_input} ${footprint_dir}/tf_list.txt
echo "TF list copied to ${footprint_dir}/tf_list.txt"
echo "TFs requested:"; cat ${footprint_dir}/tf_list.txt

# ============================================================
# EXTRACT CONSENSUS SEQUENCES FROM HOCOMOCO — generate motifs_to_pwm.tsv
# ============================================================
echo "--- Extracting consensus sequences from HOCOMOCOv12 ---"; date

conda activate fastQTL
python3 << PYEOF
import numpy as np

def pwm_to_consensus(matrix):
    bases = ['A', 'C', 'G', 'T']
    return ''.join(bases[np.argmax(row)] if max(row) >= 0.5 else 'N' for row in matrix)

# read TF list from copied file
tf_list_path = '${footprint_dir}/tf_list.txt'
with open(tf_list_path) as f:
    tf_list = [line.strip() for line in f if line.strip()]

print(f"Searching for {len(tf_list)} TFs: {tf_list}")

meme_file = '${meme_file}'

results = {}
current_motif = None
reading_matrix = False
matrix = []

with open(meme_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith('MOTIF'):
            if current_motif and matrix:
                tf_name = current_motif.split('.')[0]
                if tf_name in tf_list and tf_name not in results:
                    results[tf_name] = (current_motif, pwm_to_consensus(matrix))
            current_motif = line.split()[1]
            matrix = []
            reading_matrix = False
        elif line.startswith('letter-probability matrix'):
            reading_matrix = True
        elif reading_matrix and line and line[0].isdigit():
            matrix.append(list(map(float, line.split())))
        elif line.startswith('URL'):
            reading_matrix = False

# save last motif
if current_motif and matrix:
    tf_name = current_motif.split('.')[0]
    if tf_name in tf_list and tf_name not in results:
        results[tf_name] = (current_motif, pwm_to_consensus(matrix))

# write motifs_to_pwm.tsv and report
outfile = '${footprint_dir}/motifs_to_pwm.tsv'
found = []
not_found = []

with open(outfile, 'w') as out:
    for tf in tf_list:
        if tf in results:
            motif_name, consensus = results[tf]
            out.write(f"{tf}\t{consensus}\n")
            found.append(tf)
            print(f"  FOUND:     {tf:15s}  {consensus:30s}  ({motif_name})")
        else:
            not_found.append(tf)
            print(f"  NOT FOUND: {tf}")

print(f"\nSummary: {len(found)} found, {len(not_found)} not found")
if not_found:
    print(f"  Missing TFs: {not_found}")
    print(f"  These will be skipped — add them manually to {outfile} if needed")
print(f"\nWritten to: {outfile}")
PYEOF

# check if motifs_to_pwm.tsv was generated and has content
if [ ! -s "${footprint_dir}/motifs_to_pwm.tsv" ]; then
    echo "ERROR: motifs_to_pwm.tsv is empty or missing — check TF names match HOCOMOCOv12"
    exit 1
fi

echo "--- Motif extraction complete ---"; date
echo "motifs_to_pwm.tsv contents:"
cat ${footprint_dir}/motifs_to_pwm.tsv
conda deactivate 

# ============================================================
# WRITE INNER SINGULARITY SCRIPT
# ============================================================
rm -f ${base_dir}/run_inside_singularity_footprints.sh
cat > ${base_dir}/run_inside_singularity_footprints.sh << EOF
#!/bin/bash
base_dir=${base_dir}
data_dir=${data_dir}
footprint_dir=${footprint_dir}
prefix=${prefix}

echo "================================================"; date
echo "TF footprinting for ${prefix}"
echo "TFs to footprint:"; cat \${footprint_dir}/motifs_to_pwm.tsv | cut -f1
echo "================================================"

# ============================================================
# FOOTPRINT FOR EACH FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=\${base_dir}/fold\${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5
    fold_json=\${data_dir}/splits/fold_\${fold_n}.json

    if [ ! -f "\${model}" ]; then
        echo "WARNING: fold \${fold_n} model not found, skipping..."; continue
    fi

    fold_footprint_dir=\${footprint_dir}/fold\${fold_n}
    mkdir -p \${fold_footprint_dir}

    if [ ! -f "\${fold_footprint_dir}/footprints_done.flag" ]; then
        echo "--- Footprinting fold \${fold_n} ---"; date
        chrombpnet footprints \
            -m \${model} \
            -g \${data_dir}/downloads/hg38.fa \
            -r \${data_dir}/downloads/peaks_no_blacklist.narrowPeak \
            -fl \${fold_json} \
            -pwm_f \${footprint_dir}/motifs_to_pwm.tsv \
            -op \${fold_footprint_dir}/footprints \
            -bs 64
        touch \${fold_footprint_dir}/footprints_done.flag
        echo "--- fold \${fold_n} done ---"; date
    else
        echo "fold \${fold_n} footprints already exist, skipping..."
    fi
done

echo "All folds footprinted"; date
echo "================================================"
echo "Summary:"
echo "  Footprint files: \${footprint_dir}/fold*/footprints*"
echo "================================================"; date
EOF

# ============================================================
# RUN INSIDE SINGULARITY
# ============================================================
echo "Waiting for filesystem sync..."; date; sleep 10
singularity exec --nv ${sif} bash ${base_dir}/run_inside_singularity_footprints.sh
