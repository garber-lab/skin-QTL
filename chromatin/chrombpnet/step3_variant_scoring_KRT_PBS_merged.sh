#!/bin/bash
#BSUB -J chrombpnet_variant_score_KRT_PBS
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=15000]"
#BSUB -q gpu
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step3_variant_score_KRT_PBS.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step3_variant_score_KRT_PBS.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE
# ============================================================
BIAS_FACTOR=0.5
sample=KRT_PBS
prefix=ATAC_${sample}_merged_alldonors

base_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}/biasFactor${BIAS_FACTOR}"
data_dir="${base_dir}/data"
score_dir="${base_dir}/variant_scores"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

# paths
variant_scorer="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src"
conda_site_pkgs="/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages"

# variant file
variant_input="${score_dir}/variant_list_original_schema.txt"

# annotation files
genes_bed="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.bed"
peaks_bed="${data_dir}/downloads/peaks_no_blacklist.narrowPeak"

mkdir -p ${score_dir}


# ============================================================
# STEP 1: PREPARE VARIANT INPUT FILE FROM VCF
# ============================================================
echo "--- Preparing variant input file ---"; date

module load bcftools
module load samtools

snp_list="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/new_sigQTLs_to_score_in_chromBPNet_part3.txt"
#snp_list="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/new_sigQTLs_to_score_in_chromBPNet.tsv"
vcf="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz"
snp_list_temp="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/QTL_sig_snps_rsids.txt"
vcf_raw="${score_dir}/vcf_raw.txt"

# extract rsIDs only (column 1, no count)
awk 'NR>1' ${snp_list} | cut -f1 | sort -u > ${snp_list_temp}
echo "Total rsIDs to look up: $(wc -l < ${snp_list_temp})"

# extract matching variants from VCF — tab-separated output
#snp_list_temp="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.5/data/downloads/snp_list_for_debug.txt"
bcftools query \
    -i "ID=@${snp_list_temp}" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
    ${vcf} > ${vcf_raw}

echo "Variants found in VCF: $(wc -l < ${vcf_raw})"
echo "Preview:"; head -3 ${vcf_raw}

# ============================================================
# STEP 2: PROCESS ALLELES AND FORMAT FOR VARIANT-SCORER
# original schema: chr, pos, variant_id, allele1, allele2
# allele1 must match the hg38 reference genome
# ============================================================

echo "--- Processing alleles ---"
date

cd ${score_dir}

python3 << 'PYEOF'
import subprocess

vcf_raw = "vcf_raw.txt"
out_file = "variant_list_original_schema.txt"
meta_file = "variant_list_original_schema.metadata.tsv"
bad_file = "variant_list_original_schema.bad_ref_mismatch.tsv"

genome_fa = "/share/data/umw_biocore/genome_data/human/hg38/hg38.fa"

written = 0
skipped_mnp = 0
skipped_multiallelic = 0
skipped_ref_mismatch = 0
snp = 0
indel = 0
insertion = 0
deletion = 0
swapped = 0


def genome_seq_1based(chrom, pos, length):
    start = pos
    end = pos + length - 1
    region = f"{chrom}:{start}-{end}"

    result = subprocess.check_output(
        ["samtools", "faidx", genome_fa, region],
        text=True
    )

    seq = "".join(
        line.strip().upper()
        for line in result.splitlines()
        if not line.startswith(">")
    )

    return seq


def genome_matches(chrom, pos, allele):
    """
    Check whether allele matches hg38 at this position.
    Allele '-' cannot match genome directly.
    """
    if allele == "-":
        return False
    return genome_seq_1based(chrom, pos, len(allele)) == allele.upper()


def trim_variant(pos, ref_allele, alt_allele):
    """
    Trim shared prefix/suffix after deciding which allele is the genome reference.

    pos is kept in the same coordinate convention used by variant-scorer.
    """

    ref = ref_allele.upper()
    alt = alt_allele.upper()

    # trim shared prefix
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    # trim shared suffix
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    allele1 = ref if ref != "" else "-"
    allele2 = alt if alt != "" else "-"

    return pos, allele1, allele2


with open(vcf_raw) as f_in, \
     open(out_file, "w") as f_out, \
     open(meta_file, "w") as f_meta, \
     open(bad_file, "w") as f_bad:

    f_meta.write(
        "chrom\toriginal_pos\trsid\toriginal_ref\toriginal_alt\t"
        "ref_genome_allele\talt_nonref_allele\tswapped_against_ref\t"
        "final_pos\tvariant_id\tallele1\tallele2\tvariant_type\n"
    )

    f_bad.write(
        "chrom\tpos\trsid\tref\talt\t"
        "genome_seq_len_ref\tgenome_seq_len_alt\treason\n"
    )

    for line in f_in:
        if line.strip() == "" or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue

        chrom = fields[0]
        pos = int(fields[1])
        rsid = fields[2]
        ref = fields[3].upper()
        alt = fields[4].upper()

        # handle multi-allelic — take first ALT only
        if "," in alt:
            alt = alt.split(",")[0]
            skipped_multiallelic += 1

        # ----------------------------------------------------
        # First decide which raw allele matches hg38.
        # This is BEFORE trimming.
        # ----------------------------------------------------

        ref_matches = genome_matches(chrom, pos, ref)
        alt_matches = genome_matches(chrom, pos, alt)

        if ref_matches:
            ref_genome_allele = ref
            alt_nonref_allele = alt
            is_swapped = False

        elif alt_matches:
            ref_genome_allele = alt
            alt_nonref_allele = ref
            is_swapped = True
            swapped += 1

        else:
            skipped_ref_mismatch += 1

            seq_ref_len = genome_seq_1based(chrom, pos, len(ref))
            seq_alt_len = genome_seq_1based(chrom, pos, len(alt))

            f_bad.write(
                f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t"
                f"{seq_ref_len}\t{seq_alt_len}\t"
                f"neither_REF_nor_ALT_matches_hg38_at_raw_pos\n"
            )
            continue

        # ----------------------------------------------------
        # Now trim only after reference direction is determined.
        # allele1 should remain the genome-reference allele.
        # ----------------------------------------------------

        final_pos, allele1, allele2 = trim_variant(
            pos,
            ref_genome_allele,
            alt_nonref_allele
        )

        # ----------------------------------------------------
        # Classify final representation.
        # ----------------------------------------------------

        if allele1 != "-" and allele2 != "-" and len(allele1) == 1 and len(allele2) == 1:
            variant_type = "SNP"
            suffix = f"_{chrom}_{final_pos}_{allele1}_{allele2}"
            snp += 1

        elif allele1 != "-" and allele2 == "-":
            variant_type = "deletion"
            suffix = f"_{chrom}_{final_pos}_{allele1}_del"
            deletion += 1
            indel += 1

        elif allele1 == "-" and allele2 != "-":
            variant_type = "insertion"
            suffix = f"_{chrom}_{final_pos}_ins_{allele2}"
            insertion += 1
            indel += 1

        else:
            skipped_mnp += 1
            continue

        variant_id = rsid + suffix

        # ----------------------------------------------------
        # Write variant-scorer input: exactly 5 columns.
        # ----------------------------------------------------
        f_out.write(
            f"{chrom}\t{final_pos}\t{variant_id}\t{allele1}\t{allele2}\n"
        )

        # ----------------------------------------------------
        # Write metadata with swap status.
        # ----------------------------------------------------
        f_meta.write(
            f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t"
            f"{ref_genome_allele}\t{alt_nonref_allele}\t{is_swapped}\t"
            f"{final_pos}\t{variant_id}\t{allele1}\t{allele2}\t{variant_type}\n"
        )

        written += 1


print(f"Written: {written}")
print(f"Skipped MNP/complex: {skipped_mnp}")
print(f"Skipped reference mismatches: {skipped_ref_mismatch}")
print(f"Multi-allelic first ALT used: {skipped_multiallelic}")
print(f"SNPs: {snp}")
print(f"indels: {indel}")
print(f"insertions: {insertion}")
print(f"deletions: {deletion}")
print(f"Swapped against reference: {swapped}")
print(f"Output: {out_file}")
print(f"Metadata: {meta_file}")
print(f"Bad ref mismatch file: {bad_file}")
PYEOF


echo "Variant file preview:"
head -5 ${variant_input}

echo "Metadata preview:"
head -5 ${score_dir}/variant_list_original_schema.metadata.tsv

echo "Bad ref mismatch preview:"
head -5 ${score_dir}/variant_list_original_schema.bad_ref_mismatch.tsv

echo "Total variants:"
wc -l ${variant_input}

# ============================================================
# STEP 3: SCORE VARIANTS AND SUMMARIZE ACROSS FOLDS
# ============================================================
rm -f ${base_dir}/run_inside_singularity_variant_score.sh
cat > ${base_dir}/run_inside_singularity_variant_score.sh << EOF
#!/bin/bash
export PYTHONPATH=${conda_site_pkgs}

base_dir=${base_dir}
data_dir=${data_dir}
score_dir=${score_dir}
prefix=${prefix}
variant_input=${variant_input}
variant_scorer=${variant_scorer}
bias_model=${base_dir}/fold0/bias_model/models/${prefix}_fold0_bias.h5
peaks_bed=${peaks_bed}

echo "================================================"; date
echo "Variant scoring for ${prefix}"
echo "Total variants: \$(tail -n +2 \${variant_input} | wc -l)"
echo "================================================"

# ============================================================
# [3a] SCORE VARIANTS FOR EACH FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=\${base_dir}/fold\${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5
    fold_score_dir=\${score_dir}/fold\${fold_n}
    mkdir -p \${fold_score_dir}

    if [ ! -f "\${model}" ]; then
        echo "WARNING: fold \${fold_n} model not found, skipping..."; continue
    fi

    if [ ! -f "\${fold_score_dir}/scored_variants.variant_scores.tsv" ]; then
        echo "--- Scoring variants fold \${fold_n} ---"; date
        python \${variant_scorer}/variant_scoring.py \
            --list \${variant_input} \
            --genome \${data_dir}/downloads/hg38.fa \
            --model \${model} \
            --bias \${bias_model} \
            --out_prefix \${fold_score_dir}/scored_variants \
            --chrom_sizes \${data_dir}/downloads/hg38.chrom.sizes \
            --schema original \
            --peaks \${peaks_bed} \
            --num_shuf 10 \
            --no_hdf5 \
            --batch_size 64 \
            --random_seed 42
        echo "--- fold \${fold_n} done ---"; date
    else
        echo "fold \${fold_n} scores already exist, skipping..."
    fi
done

# ============================================================
# [3b] SUMMARIZE ACROSS FOLDS
# ============================================================
echo "--- Summarizing across folds ---"; date

score_list=""
for fold_n in 0 1 2 3 4; do
    f=\${score_dir}/fold\${fold_n}/scored_variants.variant_scores.tsv
    if [ -f "\${f}" ]; then
        score_list="\${score_list} fold\${fold_n}/scored_variants.variant_scores.tsv"
    else
        echo "WARNING: fold \${fold_n} score file not found, skipping from summary"
    fi
done

echo "Summarizing files: \${score_list}"

python \${variant_scorer}/variant_summary_across_folds.py \
    --score_dir \${score_dir} \
    --score_list \${score_list} \
    --out_prefix \${score_dir}/summary \
    --schema original

echo "--- Summary done ---"; date
echo "Summary file: \${score_dir}/summary.mean.variant_scores.tsv"
echo "Preview:"; head -3 \${score_dir}/summary.mean.variant_scores.tsv

echo "================================================"
echo "Variant scoring complete"
echo "================================================"; date
EOF

# ============================================================
# RUN INSIDE SINGULARITY
# ============================================================
echo "Waiting for filesystem sync..."; date; sleep 10
module load conda
source activate fastQTL
singularity exec --nv \
    --env PYTHONPATH=${conda_site_pkgs} \
    ${sif} bash ${base_dir}/run_inside_singularity_variant_score.sh
conda deactivate
