#!/bin/bash
# download GWAS summary statistics for each trait

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/

# Record number of cases and controls per trait
echo -e "trait\tcases\tcontrols" > ${dir}/summary_all_traits.txt
echo -e "vitiligo_verma2024\t1263\t449492" >> ${dir}/summary_all_traits.txt
echo -e "vitiligo_jin2016\t2853\t37405" >> ${dir}/summary_all_traits.txt
echo -e "psoriasis\t36466\t458078" >> ${dir}/summary_all_traits.txt
echo -e "cutaneous_lupus_erythematosus\t752\t450411" >> ${dir}/summary_all_traits.txt


# Vitiligo
# Vitiligo: 2016 Jin study (2,853 European ancestry cases and 37,405 European ancestry controls)
trait=vitiligo_jin2016
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# already downloaded before, GCST004785

# Psoriasis
# Dand 2025 largest meta analysis (36,466 European ancestry cases, 458,078 European ancestry controls)
trait=psoriasis_dand2025
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# harnomized GWAS
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90472001-GCST90473000/GCST90472771/harmonised/GCST90472771.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90472001-GCST90473000/GCST90472771/harmonised/GCST90472771.h.tsv.gz-meta.yaml

# Psoriasis vulgaris (plaque psoriasis)
# Independent of Dand 2025 — FinnGen is not a constituent cohort of GCST90472771.
# FinnGen R13, L12_PSORI_VULG (8,240 Finnish cases, 481,187 Finnish controls; N_eff = 32,405)
# Already GRCh38, tabix-indexed, beta/sebeta given directly. No liftOver, no harmonisation needed.
trait=psoriasis_vulgaris_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# FinnGen R13 summary stats + tabix index
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_PSORI_VULG.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_PSORI_VULG.gz.tbi

# atopic_dermatitis
# 2025 Nat Comm Oliva M Atopic dermatitis study (42,963 European ancestry cases, 408,472 European ancestry controls, 5,014 Asian ancestry cases, 171,135 Asian ancestry controls, 7,063 African ancestry cases, 15,879 African ancestry controls, 237 Hispanic or Latin American cases, 1,877 Hispanic or Latin American controls, 869 admixed ancestry cases, 4,917 admixed ancestry controls)
trait=atopic_dermatitis_oliva2025
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90503001-GCST90504000/GCST90503108/harmonised/GCST90503108.h.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90503001-GCST90504000/GCST90503108/harmonised/GCST90503108.h.tsv.gz-meta.yaml
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90503001-GCST90504000/GCST90503108/harmonised/GCST90503108.h.tsv.gz.tbi

# Atopic dermatitis (FinnGen study)
# FinnGen R13, L12_ATOPIC (34,368 Finnish cases, 427,815 Finnish controls; N_eff = 127,250)
# Independent cohort — FinnGen is not part of Oliva 2025 (GCST90503108/9).
# GRCh38, tabix-indexed, beta/sebeta given directly.
trait=atopic_dermatitis_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_ATOPIC.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_ATOPIC.gz.tbi

# Cutaneous lupus erythematosus
# 2024 Verma study (752 European ancestry cases, 450,411 European ancestry controls)
trait=cutaneous_lupus_erythematosus_verma2024
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476182/harmonised/GCST90476182.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476182/harmonised/GCST90476182.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476182/harmonised/GCST90476182.h.tsv.gz-meta.yaml

# Cutaneous lupus erythematosus — independent replication
# FinnGen R13, L12_LUPUS (877 Finnish cases, 462,830 controls; N_eff = 3,501)
trait=cutaneous_lupus_erythematosus_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_LUPUS.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_LUPUS.gz.tbi

# Systemic lupus erythematosus
# 2015 Nat Genet Bentham J (5,201 European ancestry cases, 9,066 European ancestry controls; N_eff = 13,220)
#   42 genome-wide-significant associations — the most of any SLE GWAS with sumstats.
#   Harmonised GRCh38 but NO .tbi (legacy filename) — index locally.
# Two Bentham findings worth noting before moving on, since both would silently corrupt results:
# - is_sorted: false — tabix will fail on it. Needs a coordinate sort before indexing (your parse_by_chr.sh splits fine, but the index step will error).
# - Only the hm_* columns are GRCh38. The plain chromosome/base_pair_location columns are the original build37 — the sibling file is literally named -build37.f.tsv.gz. Using those would put every SLE variant at hg19 coordinates against your hg38 eQTLs.
trait=systemic_lupus_erythematosus_bentham2015
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/harmonised/26502338-GCST003156-EFO_0002690.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/harmonised/26502338-GCST003156-EFO_0002690.h.tsv.gz-meta.yaml

# SLE — independent replication
# FinnGen R13, SLE_FG (1,340 Finnish cases, 368,758 controls; N_eff = 5,341)
trait=systemic_lupus_erythematosus_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_SLE_FG.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_SLE_FG.gz.tbi

# basal cell carcinoma
# 2024 Verma veteran study (31142 cases, 404406 controls, European)
trait=basal_cell_carcinoma_verma2024
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# harmonized GWAS stats, index, metadata
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475582/harmonised/GCST90475582.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475582/harmonised/GCST90475582.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475582/harmonised/GCST90475582.h.tsv.gz-meta.yaml


# squamous cell carcinoma (cutaneous)
# 2024 Science Verma A, MVP (19,217 European ancestry cases, 419,674 European ancestry controls; N_eff = 73,502)
#   PheCode 172.22 — same study/design as the basal_cell_carcinoma pick (PheCode 172.21).
trait=squamous_cell_carcinoma_verma2024
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# harmonized GWAS stats, index, metadata
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475583/harmonised/GCST90475583.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475583/harmonised/GCST90475583.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475583/harmonised/GCST90475583.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475583/harmonised/md5sum.txt

# rheumatoid arthritis
# updated
# 2022 Nat Genet study (22,350 European ancestry cases, 74,823 European ancestry controls)
trait=rheumatoid_arthritis
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132223/harmonised/GCST90132223.h.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132223/harmonised/GCST90132223.h.tsv.gz.tbi
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132223/harmonised/GCST90132223.h.tsv.gz-meta.yaml

# alopecia areata
# 2024 Science Verma A, MVP multi-ancestry (1,598 cases, 630,425 controls; N_eff = 6,376)
#   EUR 600/450,317 + AFR 627/120,788 (+ other arms).  eur_frac = 0.71
trait=alopecia_areata_verma2024
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480473/harmonised/GCST90480473.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480473/harmonised/GCST90480473.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480473/harmonised/GCST90480473.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480473/harmonised/md5sum.txt

# alopecia areata — independent replication
# FinnGen R13, L12_ALOPECAREATA (994 Finnish cases, 474,896 controls; N_eff = 3,968)
trait=alopecia_areata_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_ALOPECAREATA.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_ALOPECAREATA.gz.tbi


# melanoma
# 2024 Verma Science Veteran study 
trait=Melanomas_of_skin_dx_or_hx_verma2024
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475577/harmonised/GCST90475577.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475577/harmonised/GCST90475577.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475577/harmonised/GCST90475577.h.tsv.gz-meta.yaml


# Crohn's disease
# 2023 Nat Genet Liu Z, trans-ancestry IBD (20,873 cases, 346,719 controls; N_eff = 78,751)
#   EUR 13,501/331,263 + EAS 7,372/15,456.  eur_frac = 0.94 (European dominates both
#   cases and controls, so the trans-ancestry arm costs almost no ancestry match).
#   GRCh38, is_sorted: true, .tbi present.
trait=crohns_disease_liu2023
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446792/harmonised/GCST90446792.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446792/harmonised/GCST90446792.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446792/harmonised/GCST90446792.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446792/harmonised/md5sum.txt

# Ulcerative colitis
# 2023 Nat Genet Liu Z, trans-ancestry IBD (23,252 cases; N_eff ~80,000)
#   UC arm of the SAME study as crohns_disease (GCST90446792) — identical cohorts,
#   design and pipeline, so any Crohn's-vs-UC difference is biology, not study effect.
#   GRCh38, is_sorted: true, .tbi present.
trait=ulcerative_colitis_liu2023
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794/harmonised/GCST90446794.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794/harmonised/GCST90446794.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794/harmonised/GCST90446794.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794/harmonised/md5sum.txt


# Skin pigmentation 
# 2025 UKbiobank (415,018 European ancestry individuals)
trait=skin_pigmentation_ukbiobank2025
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691754/harmonised/GCST90691754.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691754/harmonised/GCST90691754.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691754/harmonised/GCST90691754.h.tsv.gz-meta.yaml


# Ankylosing spondylitis
# FinnGen R13, M13_ANKYLOSPON (3,875 Finnish cases, 481,187 controls; N_eff = 15,329)
#   Chosen over GCST005529 (Cortes 2013, N_eff 24,969, 47 assoc) because that study is
#   Immunochip: snpCount 128,935, targeted to ~186 immune loci with no imputation.
#   Most eGene cis-windows would contain too few variants to colocalize, so its
#   coverage gaps would masquerade as tissue specificity.
#   FinnGen is genome-wide, GRCh38, tabix-indexed, beta/sebeta direct.
trait=ankylosing_spondylitis_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_M13_ANKYLOSPON.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_M13_ANKYLOSPON.gz.tbi


# Lichen planus
# FinnGen R13, K11_LICHEN_PLANUS_WIDE (9,526 Finnish cases, 490,660 controls; N_eff = 37,378)
#   Interface dermatitis: cytotoxic T cells attack basal keratinocytes under IFN-gamma.
#   Best phenotypic match available for KRT x IFN-gamma.
#   CAVEAT: the "wide" definition POOLS oral and non-oral lichen planus. Reeve MP 2024
#   (Am J Hum Genet) shows these are genetically heterogeneous, but the non-oral-only
#   analysis (GCST90435153) was deposited WITHOUT summary statistics, so pooled is the
#   only option. Note the K11_ prefix is FinnGen's digestive-chapter tag (because of the
#   oral component), NOT an indication this isn't a skin disease.
#   The GWAS Catalog has nothing usable: every deposited LP study has 0 associations.
trait=lichen_planus_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_K11_LICHEN_PLANUS_WIDE.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_K11_LICHEN_PLANUS_WIDE.gz.tbi

# Hidradenitis suppurativa
# 2026 J Invest Dermatol Peacker BL, multi-biobank meta (3,941 cases, 1,435,603 controls; N_eff = 15,721)
#   EUR 2,864/1,315,903 + AFR 1,077/119,700.  eur_frac = 0.92
#   Follicular keratinocyte occlusion, TNF-driven (adalimumab-approved) -> KRT x TNF.
#   Chosen over FinnGen L12_HIDRADENITISSUP (N_eff 5,794) and Broadaway 2025 (N_eff ~13,244).
#   GRCh38, is_sorted: true, .tbi present.
trait=hidradenitis_suppurativa_peacker2026
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90726001-GCST90727000/GCST90726609/harmonised/GCST90726609.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90726001-GCST90727000/GCST90726609/harmonised/GCST90726609.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90726001-GCST90727000/GCST90726609/harmonised/GCST90726609.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90726001-GCST90727000/GCST90726609/harmonised/md5sum.txt

# Celiac disease
# 2025 Hum Reprod Shigesi N (7,173 European ancestry cases, 330,282 controls; N_eff = 28,082)
#   Genome-wide: 7,216,568 variants.
#   Chosen over the higher-association Immunochip studies -- Trynka 2011 (GCST005523,
#   36 assoc) and Ricano-Ponce 2019 (GCST010064) are targeted arrays with ~128-140k
#   variants concentrated in ~186 immune loci, so most eGene cis-windows would fall
#   below the <50-shared-SNP threshold and fail for array design, not biology.
#   0 curated associations here is a curation artifact: the paper is about endometriosis
#   x immune disease, so the Catalog tabulated the endometriosis results.
trait=celiac_disease_shigesi2025
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558088/harmonised/GCST90558088.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558088/harmonised/GCST90558088.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558088/harmonised/GCST90558088.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558088/harmonised/md5sum.txt


# Multiple sclerosis
# 2025 Hum Reprod Shigesi N (16,381 European ancestry cases, 343,623 controls; N_eff = 62,543)
#   Genome-wide: 7,216,598 variants. Highest N_eff of any MS study with sumstats.
#   NOT Beecham 2013 (GCST005531, 99 assoc) -- that is Immunochip, 161,311 variants.
#   The definitive IMSGC 2019 Science map (GCST009597, 8.3M variants, ~47k cases) has
#   fullPvalueSet=False -- no summary statistics deposited.
#   0 curated associations is a curation artifact (paper is about endometriosis).
#   NOTE: same source paper as celiac_disease -- not fully independent.
#   Relevance: IFN-beta is an approved MS therapy, so this maps onto the IFNB condition.
# Caveat: its causal cell types are peripheral immune cells and microglia, per the IMSGC map. NOT SKIN CELLS.
trait=multiple_sclerosis_shigesi2025
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558092/harmonised/GCST90558092.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558092/harmonised/GCST90558092.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558092/harmonised/GCST90558092.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90558001-GCST90559000/GCST90558092/harmonised/md5sum.txt

# Sjögren's syndrome
# FinnGen R13, M13_SJOGREN "Sicca syndrome [Sjögren]" (3,403 Finnish cases, 483,310 controls; N_eff = 13,517)
#   Canonical type I interferon-signature disease -> maps onto the IFNB/IFNG conditions.
#   FinnGen is 2.1x better powered than anything in the GWAS Catalog, where the best
#   option (GCST90837193, N_eff 4,901) is not even harmonised.
# Caveat: overlaps with RA/SLE symptoms. So if a locus colocalizes with Sjögren's and with RA or SLE, treat that as one finding rather than three, and don't read it as Sjögren's-specific.
#   GRCh38, tabix-indexed, beta/sebeta direct.
trait=sjogren_syndrome_finngen
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_M13_SJOGREN.gz
wget https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_M13_SJOGREN.gz.tbi

# Actinic keratosis
# 2024 Science Verma A, MVP (97,501 European ancestry cases, 314,618 controls; N_eff = 297,735)
#   PheCode 702.1. 44.3M variants, 66 genome-wide-significant associations.
#   The best-powered trait in this entire set — AK is extremely common in the VA
#   cohort (older, sun-exposed), so the PheCode definition works here where it fails
#   for rare diseases like vitiligo and alopecia.
#   Keratinocyte precancer: bridges normal keratinocyte biology and squamous_cell_carcinoma
#   (same MVP study family, PheCode 172.22).
#   GRCh38, is_sorted: true, .tbi present.
trait=actinic_keratosis_verma2024
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476204/harmonised/GCST90476204.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476204/harmonised/GCST90476204.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476204/harmonised/GCST90476204.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476204/harmonised/md5sum.txt

# Asthma. share FLG/barrier and Th2 genetics with Atopic Dermatitis
# autoimmune thyroid disease. strongest genetic correlate of vitiligo
# bullous pemphigoid / pemphigus. autoimmune blistering diseases
# T1D. Strong IFN component
# Urticaria. Mast-cell/skin
# Dupuytren's contracture. fibroblast proliferation with no immune component, paired against systemic sclerosis
# ****************************************************************************************
# Systemic sclerosis
# 2019 Nat Commun López-Isac E, meta-analysis of 14 GWAS (9,095 European cases, 17,584 controls; N_eff = 23,978)
#   Best available by far — every alternative is a PheCode study with <500 cases and 0 associations.
#   NOT harmonised. Bespoke parse required:
#     - OR only, no SE  -> beta = log(OR); se = |beta| / |qnorm(P/2)|
#     - A2 is '?' for some variants -> recover from reference or drop
#     - hg19 (rs3094315 at chr1:752566) -> liftOver to hg38
#     - the N column is the number of contributing cohorts, NOT sample size
trait=systemic_sclerosis_lopez2019
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009131/Lopez-Isac_prePMID_META_GWAS_SSc.meta.txt
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009131/Lopez-Isac_Readme.txt
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009131/md5sum.txt



