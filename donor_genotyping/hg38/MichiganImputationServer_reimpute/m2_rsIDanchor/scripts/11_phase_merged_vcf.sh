#!/bin/bash
# shuo.shan@umassmed.edu
#BSUB -J phase_merged_vcf[1-22]
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
#BSUB -o step11_phaseVCF_%I.out
#BSUB -e step11_phaseVCF_%I.err
#BSUB -q short
#BSUB -W 4:00
#BSUB -n 4

# Phase the merged + filtered imputed genotypes (step 9 output, per chromosome) with Beagle 5,
# using the phased 1000G 30x panel as reference. One job per autosome (chr1-22).
#
# 1. fix REF/ALT against hg38 (bcftools +fixref -m flip): ~25% of SNPs have REF/ALT swapped
#    relative to hg38 after merging the three MIS panels; genotypes are recoded accordingly
# 2. keep biallelic sites whose REF matches hg38
# 3. subset the 1000G panel to our sites (phasing only needs those -> faster, less memory)
# 4. Beagle: phase only (impute=false); sites not in the 1000G panel are dropped by Beagle
# 5. index + summary line in log/step11_summary.txt
#
# NOTE: after step 1 the allele coding follows hg38 (REF = hg38 base). For the swapped SNPs this is
# the OPPOSITE of the coding used for eQTL mapping, so per-ALT-allele effect signs differ for those SNPs.

# Load necessary modules
module load bcftools
module load htslib
module load openjdk/20.0.1
if ! command -v java > /dev/null; then echo "ERROR java not found: load your cluster's java module in this script"; exit 1; fi
java -version 2>&1 | head -1

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor
VCF_DIR=${DIR}/MIS_results/combined_output
OUT_DIR=${VCF_DIR}/final_output/phased
TMP_DIR=${OUT_DIR}/tmp
FASTA=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_FRB_PBS_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa
REF_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/literature/1000G/GRCh38_30x_phased/vcf
BEAGLE_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/literature/beagle
BEAGLE_JAR=${BEAGLE_DIR}/beagle.27Feb25.75f.jar
JAVA_MEM=28g                              # a bit below the job memory (4 cores x 8000 MB)

CHROM_LIST=({1..22})
chr_name="chr${CHROM_LIST[$((LSB_JOBINDEX-1))]}"
mkdir -p ${OUT_DIR} ${TMP_DIR} ${DIR}/log
cd ${OUT_DIR}

VCF2=${VCF_DIR}/${chr_name}_imputed_combined_hg38_rsID_filtered.vcf.gz      # step 9 output (per chromosome)
VCF_FIX=${TMP_DIR}/${chr_name}_fixref.vcf.gz
VCF_REF=${TMP_DIR}/${chr_name}_1000G_at_our_sites.vcf.gz
MAP=${BEAGLE_DIR}/chr_in_chrom_field/plink.chr${chr_name}.GRCh38.map     # Beagle map files are named plink.chrchrN...
MAP_USE=${TMP_DIR}/${chr_name}.map
PANEL=${REF_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr_name}.filtered.shapeit2-duohmm-phased.vcf.gz
OUT_PREFIX=${OUT_DIR}/${chr_name}_imputed_combined_hg38_rsID_filtered_phased

# remove intermediate files however the job ends (success, error, or LSF kill at the wall-time limit);
# the small logs (tmp/chrN_fixref.log, the Beagle .log) are kept
cleanup() {
    rm -f ${VCF_FIX} ${VCF_FIX}.tbi ${VCF_REF} ${VCF_REF}.tbi ${TMP_DIR}/${chr_name}_sites.txt ${MAP_USE}
}
trap cleanup EXIT
trap 'exit 1' TERM INT

for f in ${VCF2} ${FASTA} ${PANEL} ${MAP} ${BEAGLE_JAR}; do
    if [ ! -e ${f} ]; then echo "ERROR missing file: ${f}"; exit 1; fi
done

# ----------- 1-2. fix REF/ALT against hg38, keep biallelic REF-matching sites --------------
echo "Fixing REF/ALT for ${chr_name}..."; date
bcftools +fixref ${VCF2} -Ou -- -f ${FASTA} -m flip -d 2> ${TMP_DIR}/${chr_name}_fixref.log |\
        bcftools view -m2 -M2 -Ou |\
        bcftools norm -f ${FASTA} -c x -Oz -o ${VCF_FIX}
tabix -f ${VCF_FIX}

# ----------- 3. 1000G panel at our sites only --------------
echo "Subsetting 1000G panel..."; date
bcftools query -f '%CHROM\t%POS\n' ${VCF_FIX} > ${TMP_DIR}/${chr_name}_sites.txt
bcftools view -T ${TMP_DIR}/${chr_name}_sites.txt -m2 -M2 ${PANEL} -Oz -o ${VCF_REF}
tabix -f ${VCF_REF}

# genetic map: Beagle needs the same chromosome names as the VCF (chrN)
awk -v c=${chr_name} 'BEGIN{OFS="\t"} {$1 = c; print}' ${MAP} > ${MAP_USE}

# ----------- 4. phase with Beagle --------------
echo "Phasing ${chr_name} with Beagle..."; date
java -Xmx${JAVA_MEM} -jar ${BEAGLE_JAR} \
        gt=${VCF_FIX} ref=${VCF_REF} map=${MAP_USE} chrom=${chr_name} \
        impute=false nthreads=${LSB_DJOB_NUMPROC:-4} out=${OUT_PREFIX}
if [ ! -s ${OUT_PREFIX}.vcf.gz ]; then echo "ERROR Beagle failed, see ${OUT_PREFIX}.log"; exit 1; fi
tabix -f ${OUT_PREFIX}.vcf.gz

# ----------- 5. summary --------------
n1=$(bcftools view -H ${VCF2} | wc -l | cut -d' ' -f1)
n2=$(bcftools view -H ${VCF_FIX} | wc -l | cut -d' ' -f1)
n3=$(bcftools view -H ${OUT_PREFIX}.vcf.gz | wc -l | cut -d' ' -f1)
n_unphased=$(bcftools query -f '[%GT\n]' ${OUT_PREFIX}.vcf.gz | grep -c '/')
pct=$(echo "${n3}/${n1}*100" | bc -l)
pct=$(printf "%.2f" ${pct})

# chr  input  after_fixref  phased  pct_kept  unphased_GTs
echo -e "${chr_name}\t${n1}\t${n2}\t${n3}\t${pct}\t${n_unphased}" >> ${DIR}/log/step11_summary.txt

echo "Done ${chr_name}: ${OUT_PREFIX}.vcf.gz"; date

