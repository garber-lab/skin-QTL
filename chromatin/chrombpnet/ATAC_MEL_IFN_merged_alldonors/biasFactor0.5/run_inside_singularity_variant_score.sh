#!/bin/bash
export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages

base_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5
data_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/data
score_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/variant_scores
prefix=ATAC_MEL_IFN_merged_alldonors
variant_input=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/variant_scores/variant_list_original_schema.txt
variant_scorer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src
bias_model=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/fold0/bias_model/models/ATAC_MEL_IFN_merged_alldonors_fold0_bias.h5
peaks_bed=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/data/downloads/peaks_no_blacklist.narrowPeak

echo "================================================"; date
echo "Variant scoring for ATAC_MEL_IFN_merged_alldonors"
echo "Total variants: $(tail -n +2 ${variant_input} | wc -l)"
echo "================================================"

# ============================================================
# [3a] SCORE VARIANTS FOR EACH FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=${base_dir}/fold${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5
    fold_score_dir=${score_dir}/fold${fold_n}
    mkdir -p ${fold_score_dir}

    if [ ! -f "${model}" ]; then
        echo "WARNING: fold ${fold_n} model not found, skipping..."; continue
    fi

    if [ ! -f "${fold_score_dir}/scored_variants.variant_scores.tsv" ]; then
        echo "--- Scoring variants fold ${fold_n} ---"; date
        python ${variant_scorer}/variant_scoring.py             --list ${variant_input}             --genome ${data_dir}/downloads/hg38.fa             --model ${model}             --bias ${bias_model}             --out_prefix ${fold_score_dir}/scored_variants             --chrom_sizes ${data_dir}/downloads/hg38.chrom.sizes             --schema original             --peaks ${peaks_bed}             --num_shuf 10             --no_hdf5             --batch_size 64             --random_seed 42

        echo "--- fold ${fold_n} done ---"; date
    else
        echo "fold ${fold_n} scores already exist, skipping..."
    fi
done

# ============================================================
# [3b] SUMMARIZE ACROSS FOLDS
# ============================================================
echo "--- Summarizing across folds ---"; date

score_list=""
for fold_n in 0 1 2 3 4; do
    f=${score_dir}/fold${fold_n}/scored_variants.variant_scores.tsv
    if [ -f "${f}" ]; then
        score_list="${score_list} fold${fold_n}/scored_variants.variant_scores.tsv"
    else
        echo "WARNING: fold ${fold_n} score file not found, skipping from summary"
    fi
done

echo "Summarizing files: ${score_list}"

python ${variant_scorer}/variant_summary_across_folds.py     --score_dir ${score_dir}     --score_list ${score_list}     --out_prefix ${score_dir}/summary     --schema original

echo "--- Summary done ---"; date
echo "Summary file: ${score_dir}/summary.mean.variant_scores.tsv"
echo "Preview:"; head -3 ${score_dir}/summary.mean.variant_scores.tsv

echo "================================================"
echo "Variant scoring complete"
echo "================================================"; date
