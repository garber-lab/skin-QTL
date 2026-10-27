#!/bin/bash
# shuo.shan@umassmed.edu
# Pre-processing for the MATCHED-DONOR (IFNG & IFNB) QTL pipeline.
#
# Gene universe = UNION of the expressed-&-induced gene lists for IFNG and IFNB,
# across all 3 celltypes (i.e. the genes we actually care about for the IFNG-vs-IFNB
# matched-donor comparison). We build our OWN chunks over that universe at
# 25 genes/chunk: gene chunks, variant-gene pairs per chunk, genotype per chunk,
# and the VST subset per chunk. Everything lands under this pipeline's data/chunks.
#
# Reused from the genome-wide run (read-only): the gene-TSS bed and the per-context
# expressed/induced gene lists (gene_list_<ct>_<cond>_induced.txt). Nothing here
# writes into the genome-wide tree.
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J subsetPrep
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step0_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step0_%J_%I.err"

set -e

# ----------------------------------------- #
# -------- Global Variables        -------- #
# ----------------------------------------- #
GW_ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models              # genome-wide run: reuse gene lists + gene TSS
ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB   # this pipeline's outputs
SCRIPTDIR=${ROOT}/scripts
DATADIR=${ROOT}/data
CHUNK_DIR=${DATADIR}/chunks
META_FILE=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt
VST_FILE=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz

genemodel=${GW_ROOT}/data/gencode_v34_gene_TSS.bed                          # reused (bed of gene TSS)
genelist=${DATADIR}/union_IFNG_IFNB_induced_genes.txt                       # this run's gene universe
gene_tss_subset=${DATADIR}/union_gene_TSS.bed
pairs=${DATADIR}/SNPs_near_TSS.bed                                          # union-only variant-gene pairs (master)

GENES_PER_CHUNK=25

# ----------------------------------------- #
# -------- Make output tree        -------- #
# ----------------------------------------- #
mkdir -p ${DATADIR} ${CHUNK_DIR} ${ROOT}/scripts ${ROOT}/logs ${ROOT}/results
for ct in FRB KRT MEL; do
	mkdir -p ${ROOT}/${ct}/results ${ROOT}/${ct}/results_QC ${ROOT}/${ct}/eigenMT/results
done

# ----------------------------------------- #
# -------- Genotyped donor + subset ------- #
# ----------------------------------------- #
module load bcftools
module load bedtools
module load htslib

# Clean sample names identically to the genotype-chunk step so ids share meta$donor's namespace.
echo "Building genotyped donor list from ${vcf}"; date
bcftools query -l ${vcf} | cut -d'_' -f1 | sed 's/skineQTL-//' | sed 's/F0/F/' | sort -u \
	> ${DATADIR}/genotyped_donors.txt
echo "Genotyped donors: $(wc -l < ${DATADIR}/genotyped_donors.txt)"

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step0.0_define_donor_subset.R \
		${ROOT} ${META_FILE} ${DATADIR}/genotyped_donors.txt

# ----------------------------------------- #
# -------- Build union gene set   --------- #
# ----------------------------------------- #
# Union of IFNG-induced and IFNB-induced (= expressed & induced) gene lists over
# all 3 celltypes. These lists come from the genome-wide step0.1.
echo "Building union of IFNG & IFNB induced genes across celltypes..."; date
: > ${genelist}.tmp
for ct in FRB KRT MEL; do
	for cond in IFNG IFNB; do
		src=${GW_ROOT}/data/gene_list_${ct}_${cond}_induced.txt
		if [ ! -s "${src}" ]; then
			echo "ERROR: missing induced gene list ${src} -- run the genome-wide step0.1 first." >&2
			exit 1
		fi
		cat "${src}" >> ${genelist}.tmp
	done
done
sort -u ${genelist}.tmp > ${genelist}
rm ${genelist}.tmp
echo "Union gene universe: $(wc -l < ${genelist}) genes"

# ----------------------------------------- #
# -------- Make SNP-gene pair list -------- #
# ----------------------------------------- #
# Restrict the gene-TSS bed to the union genes, then window +/-500kb against all
# genotyped SNPs. Filtering -a to union genes first keeps bedtools fast.
awk 'NR==FNR{g[$1]; next} ($4 in g)' ${genelist} ${genemodel} \
	| sort -k1,1 -k2,2n > ${gene_tss_subset}

# all genotyped SNPs as bed (reuse genome-wide file if present, else build)
if [ -s "${GW_ROOT}/data/all_genotyped_SNPs.bed" ]; then
	all_snps=${GW_ROOT}/data/all_genotyped_SNPs.bed
	echo "Reusing genome-wide all_genotyped_SNPs.bed"
else
	all_snps=${DATADIR}/all_genotyped_SNPs.bed
	echo "Building all_genotyped_SNPs.bed from VCF"; date
	bcftools view -H ${vcf} | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ${all_snps}
fi

# header + windowed pairs (11 cols: gene TSS 7 + SNP bed 4)
echo -e "gene_chr\tgene_start\tgene_end\tgene_name\t0\tstrand\tgene_id\tSNP_chr\tSNP_start\tSNP_end\tSNP_ID" > ${pairs}
bedtools window -w 500000 -a ${gene_tss_subset} -b ${all_snps} >> ${pairs}
echo "Variant-gene pairs (union): $(($(wc -l < ${pairs}) - 1))"; date

# ----------------------------------------- #
# ------ chunk genes (25 / chunk) --------- #
# ----------------------------------------- #
echo "Chunking union genes at ${GENES_PER_CHUNK} genes/chunk..."; date
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step0.2_compile_gene_chunks.R ${pairs} ${CHUNK_DIR} ${GENES_PER_CHUNK}
NCHUNK=$(ls -1 ${CHUNK_DIR}/gene_chunk_* | wc -l)
echo "Num gene chunks: ${NCHUNK}"

# ----- for each gene chunk, extract matching pair rows (with header) ----- #
header=$(head -n1 ${pairs})
for gfile in ${CHUNK_DIR}/gene_chunk_*; do
	chunk_id=$(basename ${gfile} | sed 's/gene_chunk_//')
	out_pairs=${CHUNK_DIR}/pairs_chunk_${chunk_id}.tsv
	echo "Building ${out_pairs}"; date
	awk -F'\t' -v OFS='\t' -v header="${header}" '
		NR==FNR { genes[$1]=1; next }
		FNR==1 { print header; next }
		($4 in genes) { print }
	' "${gfile}" "${pairs}" > "${out_pairs}"
done
echo "Done gene-chunked pair files in ${CHUNK_DIR}/"

# ----------------------------------------- #
# ------ chunk VST per gene chunk --------- #
# ----------------------------------------- #
echo "Chunking VST to match gene chunks..."; date
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step0.4_chunk_vst.R ${VST_FILE} ${CHUNK_DIR}

# -------------------------------------------------- #
# ----- Compile Genotype table for each chunk ------ #
# -------------------------------------------------- #
# Genotype tables keep ALL genotyped donors (cleaned names); the matched-donor
# restriction is applied later at model-fit time (step1), not here.
for chunk in ${CHUNK_DIR}/pairs_chunk_*.tsv; do
	chunk_id=$(basename ${chunk} | sed 's/pairs_chunk_//' | sed 's/\.tsv//g')
	chunk_base=pairs_chunk_${chunk_id}

	snplist=${CHUNK_DIR}/temp.${chunk_base}.snps.txt
	vcf_subset=${CHUNK_DIR}/temp.${chunk_base}.vcf.gz
	samples=${CHUNK_DIR}/temp.${chunk_base}.samples.txt
	header_f=${CHUNK_DIR}/temp.${chunk_base}.header.txt
	body=${CHUNK_DIR}/temp.${chunk_base}.body.tsv
	geno=${CHUNK_DIR}/temp.${chunk_base}.genotype.tsv
	geno_num=${CHUNK_DIR}/genotype_${chunk_base}.tsv

	echo "Subsetting genotype VCF for ${chunk}"; date
	cd ${CHUNK_DIR}

	# 1) collect SNP IDs from this chunk
	awk 'NR>1 {print $11}' "$chunk" | sort -u > "$snplist"

	# 2) extract VCF records
	bcftools view --threads 8 --include ID==@"$snplist" "$vcf" -Oz -o "$vcf_subset"

	# 3) clean sample names
	bcftools query -l ${vcf_subset} | cut -d'_' -f1 | sed 's/skineQTL-//' | sed 's/F0/F/' > ${samples}

	# 4) build header
	printf "CHROM\tPOS\tID\tREF\tALT\t" > "$header_f"
	paste -sd '\t' "$samples" >> "$header_f"

	# 5) extract GT matrix
	bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf_subset" > "$body"
	cat "$header_f" "$body" > "$geno"

	# 6) convert genotypes to numeric
	sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' "$geno" > "$geno_num"

	# 7) cleanup intermediates
	rm "$snplist" "$vcf_subset" "$samples" "$header_f" "$body" "$geno"
done

echo "Done. Chunks (gene / pairs / genotype / vst) in ${CHUNK_DIR}/"; date
echo ""
echo ">>> IMPORTANT: set the job-array upper bound to the chunk count in step1/2/3 submit scripts:"
echo "    #BSUB -J QTL[1-${NCHUNK}]%70   (and QC[1-${NCHUNK}], eigenMT[1-${NCHUNK}])"
